//! Simple IR for atom/bond-based molecular models.
//!
//! This module provides two molecule representations:
//! - `Molecule`: Basic molecule with atoms, bonds, and rings
//! - `ExtendedMolecule`: Full molecule with MDL extensions (SGroups, RGroups, etc.)

use std::collections::BTreeMap;
use std::fmt::{self, Display};
use std::mem;

use serde::{Deserialize, Serialize};
use umol_data::{Element, NamedIsotope};

use crate::position::Point3D;
use crate::span::Span;

/// Input molecular format
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SourceFormat {
    MOL,
    SMILES,
    SMARTS,
    #[default]
    UNKNOWN,
}

/// Basic molecule IR - core graph data only
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub rings: Vec<Ring>,
    pub source_format: SourceFormat,
}

impl Molecule {
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Get sum formula in Hill notation (C first, H second, then alphabetically)
    pub fn sum_formula(&self) -> String {
        use std::collections::BTreeMap;
        use umol_data::Element;

        let mut atom_counts: BTreeMap<[u8; 2], (Element, usize)> = BTreeMap::new();
        let mut c_count = 0usize;
        let mut h_count = 0usize;
        let mut charge = 0i32;

        for atom in &self.atoms {
            let element = match &atom.symbol {
                AtomSymbol::Element(e) => Some(*e),
                AtomSymbol::NamedIsotope(i) => Some(i.element()),
                _ => None,
            };
            if let Some(element) = element {
                match element {
                    e if e == Element::C => c_count += 1,
                    e if e == Element::H => h_count += 1,
                    e => {
                        let key = element_symbol_key(e);
                        atom_counts.entry(key).or_insert((e, 0)).1 += 1;
                    }
                }
            }
            if let Some(ch) = atom.charge {
                charge += ch as i32;
            }
        }

        format_sum_formula(c_count, h_count, atom_counts, charge)
    }
}

/// Extended molecule IR - includes MDL extensions (SGroups, RGroups, etc.)
/// This is a flat structure using ExtendedAtom and ExtendedBond.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ExtendedMolecule {
    // Core structure with extended atom/bond types
    pub atoms: Vec<ExtendedAtom>,
    pub bonds: Vec<ExtendedBond>,
    pub rings: Vec<Ring>,
    pub source_format: SourceFormat,

    // MDL extensions
    pub sgroups: BTreeMap<usize, SGroup>,
    pub rgroups: BTreeMap<usize, RGroup>,

    // Additional structure
    pub fragments: Vec<Fragment>,
    pub links: Vec<Link>,
    pub electrons: Option<u32>,

    // Properties and metadata
    pub properties: Vec<Property>,
    pub comments: Vec<String>,
}

impl ExtendedMolecule {
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Convert to basic Molecule (drops CTAB-specific fields, S-groups, R-groups)
    pub fn to_basic(&self) -> Molecule {
        Molecule {
            atoms: self.atoms.iter().map(|a| a.to_basic()).collect(),
            bonds: self.bonds.iter().map(|b| b.to_basic()).collect(),
            rings: self.rings.clone(),
            source_format: self.source_format,
        }
    }

    /// Get sum formula in Hill notation (C first, H second, then alphabetically)
    pub fn sum_formula(&self) -> String {
        use std::collections::BTreeMap;
        use umol_data::Element;

        let mut atom_counts: BTreeMap<[u8; 2], (Element, usize)> = BTreeMap::new();
        let mut c_count = 0usize;
        let mut h_count = 0usize;
        let mut charge = 0i32;

        for atom in &self.atoms {
            let element = match &atom.symbol {
                AtomSymbol::Element(e) => Some(*e),
                AtomSymbol::NamedIsotope(i) => Some(i.element()),
                _ => None,
            };
            if let Some(element) = element {
                match element {
                    e if e == Element::C => c_count += 1,
                    e if e == Element::H => h_count += 1,
                    e => {
                        let key = element_symbol_key(e);
                        atom_counts.entry(key).or_insert((e, 0)).1 += 1;
                    }
                }
            }
            if let Some(ch) = atom.charge {
                charge += ch as i32;
            }
        }

        format_sum_formula(c_count, h_count, atom_counts, charge)
    }

    /// Extract basic molecule (converts ExtendedAtom/ExtendedBond to basic types)
    pub fn to_molecule(&self) -> Molecule {
        Molecule {
            atoms: self.atoms.iter().map(|a| a.to_basic()).collect(),
            bonds: self.bonds.iter().map(|b| b.to_basic()).collect(),
            rings: self.rings.clone(),
            source_format: self.source_format,
        }
    }

    /// Create from basic molecule
    pub fn from_molecule(mol: Molecule) -> Self {
        Self {
            atoms: mol.atoms.into_iter().map(ExtendedAtom::from_basic).collect(),
            bonds: mol.bonds.into_iter().map(ExtendedBond::from_basic).collect(),
            rings: mol.rings,
            source_format: mol.source_format,
            ..Default::default()
        }
    }
}

impl From<Molecule> for ExtendedMolecule {
    fn from(mol: Molecule) -> Self {
        ExtendedMolecule::from_molecule(mol)
    }
}
/// Atom IR
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    pub symbol: AtomSymbol,
    pub position: Option<Point3D>,
    pub charge: Option<i8>,
    pub isotope: Option<u32>,
    pub radical: Option<AtomRadical>,
    pub hydrogens: Option<u8>,
    pub implicit_h: bool,
    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,
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
    AtomList(AtomList),
    Variable(Variable),
    RGroup(RGroup),
    LonePair,
    Pseudoatom(String),
    #[default]
    Unknown,
}

impl AtomSymbol {
    /// Returns true if this is an extended atom structure (not a simple element or isotope)
    pub fn is_extended(&self) -> bool {
        !matches!(self, AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_))
    }
}

/// Variable atom
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Variable {}

/// Extended query atom types (superset of MOL + SMARTS)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
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

impl QueryAtom {
    pub fn symbol(&self) -> &str {
        match self {
            QueryAtom::Any => "*",
            QueryAtom::Heavy => "A",
            QueryAtom::Heteroatom => "Q",
            QueryAtom::Halogen => "X",
            QueryAtom::Metal => "M",
            QueryAtom::HeavyOrH => "AH",
            QueryAtom::HeteroatomOrH => "QH",
            QueryAtom::HalogenOrH => "XH",
            QueryAtom::MetalOrH => "MH",
            QueryAtom::Unknown => "?",
        }
    }

    pub fn from_symbol_bytes(s: &[u8]) -> Option<QueryAtom> {
        match s {
            b"*" => Some(QueryAtom::Any),
            b"A" => Some(QueryAtom::Heavy),
            b"Q" => Some(QueryAtom::Heteroatom),
            b"X" => Some(QueryAtom::Halogen),
            b"M" => Some(QueryAtom::Metal),
            b"AH" => Some(QueryAtom::HeavyOrH),
            b"QH" => Some(QueryAtom::HeteroatomOrH),
            b"XH" => Some(QueryAtom::HalogenOrH),
            b"MH" => Some(QueryAtom::MetalOrH),
            _ => None,
        }
    }

    pub fn from_symbol_str(s: &str) -> Option<QueryAtom> {
        Self::from_symbol_bytes(s.as_bytes())
    }
}

impl Display for QueryAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

/// Atom list (inclusion or exclusion list of elements)
#[derive(Debug, Default, Clone, PartialEq, Serialize, Deserialize)]
pub struct AtomList {
    pub elements: Vec<Element>,
    pub exclusion: bool,
}

/// Radical type (unified from CTAB and SMILES)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomRadical {
    Singlet,   // 0 unpaired electrons, but still a radical center
    Doublet,   // 1 unpaired electron
    Triplet,   // 2 unpaired electrons
    Other(u8), // Other radical states (for extensibility)
}

impl AtomRadical {
    /// Convert from unpaired electron count
    pub fn from_unpaired_e(n: u8) -> Option<Self> {
        match n {
            0 => None,
            1 => Some(AtomRadical::Doublet),
            2 => Some(AtomRadical::Triplet),
            n => Some(AtomRadical::Other(n)),
        }
    }

    /// Convert to unpaired electron count
    pub fn to_unpaired_e(self) -> u8 {
        match self {
            AtomRadical::Singlet => 0,
            AtomRadical::Doublet => 1,
            AtomRadical::Triplet => 2,
            AtomRadical::Other(n) => n,
        }
    }
}

/// Chirality
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Chirality {
    Clockwise,
    CounterClockwise,
    Tetrahedral { arr: u32 },
    Allenal { arr: u32 },
    SquarePlanar { arr: u32 },
    TrigonalBipyramidal { arr: u32 },
    Octahedral { arr: u32 },
    #[default]
    Unknown,
}

/// Extended atom IR - includes all CTAB-specific fields (flat structure)
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ExtendedAtom {
    // Core fields (same as Atom)
    pub symbol: AtomSymbol,
    pub position: Option<Point3D>,
    pub charge: Option<i8>,
    pub isotope: Option<u32>,
    pub radical: Option<AtomRadical>,
    pub hydrogens: Option<u8>,
    pub implicit_h: bool,
    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,
    pub span: Option<Span>,

    // CTAB-specific fields
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
    pub properties: std::collections::HashMap<String, String>,
}

impl ExtendedAtom {
    pub fn new(symbol: AtomSymbol) -> Self {
        Self {
            symbol,
            ..Default::default()
        }
    }

    /// Convert to basic Atom (drops CTAB-specific fields)
    pub fn to_basic(&self) -> Atom {
        Atom {
            symbol: self.symbol.clone(),
            position: self.position.clone(),
            charge: self.charge,
            isotope: self.isotope,
            radical: self.radical,
            hydrogens: self.hydrogens,
            implicit_h: self.implicit_h,
            aromatic: self.aromatic,
            chirality: self.chirality,
            class: self.class,
            span: self.span,
        }
    }

    /// Create from basic Atom
    pub fn from_basic(atom: Atom) -> Self {
        Self {
            symbol: atom.symbol,
            position: atom.position,
            charge: atom.charge,
            isotope: atom.isotope,
            radical: atom.radical,
            hydrogens: atom.hydrogens,
            implicit_h: atom.implicit_h,
            aromatic: atom.aromatic,
            chirality: atom.chirality,
            class: atom.class,
            span: atom.span,
            ..Default::default()
        }
    }
}

/// Atom stereo parity (from CTAB)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomStereoParity {
    Odd,   // Clockwise / R
    Even,  // Counter-Clockwise / S
    Either,
}

/// Atom stereo care (from CTAB)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomStereoCare {
    Care, // Stereo should be considered
}

/// Atom inversion/retention (from CTAB reactions)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomInversionRetention {
    Inverted,
    Retained,
}

/// Atom exact change flag (from CTAB reactions)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomExactChange {
    Match,
}

/// Attachment point type (for R-groups)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AttachmentPointType {
    First,
    Second,
    Both,
}

/// Ring bond count query (from CTAB)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RingBondCount {
    AsDrawn,
    NoRingBonds,
    R2,
    R3,
    R4Plus,
}

/// Substitution count query (from CTAB)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
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
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct UnsaturatedAtom;

/// Link atom (for polymers)
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct LinkAtom {
    pub repeat_count: u8,
    pub subs_index1: usize,
    pub subs_index2: Option<usize>,
}

/// Bond IR
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    pub start_atom: u32,
    pub end_atom: u32,
    pub order: BondOrder,
    pub ring: Option<u32>,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDir>,
    pub span: Option<Span>,
}

impl Bond {
    pub fn new(start_atom: u32, end_atom: u32, order: BondOrder) -> Self {
        Self {
            start_atom,
            end_atom,
            order,
            ..Default::default()
        }
    }

    pub fn from_order(start_atom: u32, end_atom: u32, order: BondOrder) -> Self {
        Self::new(start_atom, end_atom, order)
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
            order,
            span: Span::from_bytes_opt(span_start, span_end),
            ..Default::default()
        }
    }
}

/// Unified bond order (merges simple_ir::BondOrder and ctab::BondType)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondOrder {
    // Standard bond orders
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Sextuple,
    Aromatic,
    // Query bond types
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
    Any,
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
            BondOrder::Quintuple => "~5",
            BondOrder::Sextuple => "~6",
            BondOrder::Aromatic => ":",
            BondOrder::SingleOrDouble => "~",
            BondOrder::SingleOrAromatic => "~",
            BondOrder::DoubleOrAromatic => "~",
            BondOrder::Any => "~",
            BondOrder::Unknown => "?",
        }
    }

    pub fn is_query(&self) -> bool {
        matches!(
            self,
            BondOrder::SingleOrDouble
                | BondOrder::SingleOrAromatic
                | BondOrder::DoubleOrAromatic
                | BondOrder::Any
        )
    }

    pub fn is_extended(&self) -> bool {
        matches!(
            self,
            BondOrder::Zero | BondOrder::Quadruple | BondOrder::Quintuple | BondOrder::Sextuple
        )
    }
}

/// Single bond direction/wedging
/// In MOL files: Up=Wedge (code 1), Down=Dash (code 6)
/// In SMILES: Up=/, Down=\
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondDir {
    Up,    // MOL: Wedge (code 1), SMILES: /
    Down,  // MOL: Dash (code 6), SMILES: \
    #[default]
    Either, // MOL code 4 (Either)
}

impl BondDir {
    pub fn is_default(&self) -> bool {
        matches!(self, d if *d == Default::default())
    }
}

/// Double-bond stereochemistry (E/Z) annotation in IR
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondStereo {
    Cis,
    Trans,
    #[default]
    Either,
}

/// Extended bond IR - includes all CTAB-specific fields (flat structure)
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ExtendedBond {
    // Core fields (same as Bond)
    pub start_atom: u32,
    pub end_atom: u32,
    pub order: BondOrder,
    pub ring: Option<u32>,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDir>,
    pub span: Option<Span>,

    // CTAB-specific fields
    pub topology: Option<BondTopology>,
    pub reacting_center: Option<BondReactingCenter>,
    pub properties: std::collections::HashMap<String, String>,
}

impl ExtendedBond {
    pub fn new(start_atom: u32, end_atom: u32, order: BondOrder) -> Self {
        Self {
            start_atom,
            end_atom,
            order,
            ..Default::default()
        }
    }

    /// Create a bond with just the order (atom indices set to 0, for later update)
    pub fn with_order(order: BondOrder) -> Self {
        Self {
            start_atom: 0,
            end_atom: 0,
            order,
            ..Default::default()
        }
    }

    /// Convert to basic Bond (drops CTAB-specific fields)
    pub fn to_basic(&self) -> Bond {
        Bond {
            start_atom: self.start_atom,
            end_atom: self.end_atom,
            order: self.order,
            ring: self.ring,
            stereo: self.stereo,
            direction: self.direction,
            span: self.span,
        }
    }

    /// Create from basic Bond
    pub fn from_basic(bond: Bond) -> Self {
        Self {
            start_atom: bond.start_atom,
            end_atom: bond.end_atom,
            order: bond.order,
            ring: bond.ring,
            stereo: bond.stereo,
            direction: bond.direction,
            span: bond.span,
            ..Default::default()
        }
    }
}

/// Bond topology (chain, ring, either) query
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondTopology {
    Chain, // MOL code 2
    Ring,  // MOL code 1
    #[default]
    Either, // MOL code 0 (default/unspecified)
}

impl BondTopology {
    pub fn is_default(&self) -> bool {
        matches!(self, t if *t == Default::default())
    }
}

bitflags::bitflags! {
    /// Bond reacting center (from CTAB reactions) - bitflags
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
    pub struct BondReactingCenter: i16 {
        const UNMARKED         = 0b00000000;
        const CENTER           = 0b00000001;
        const NOT_CENTER       = 0b00000010;
        const NO_CHANGE        = 0b00000100;
        const MADE_BROKEN      = 0b00001000;
        const ORDER_CHANGED    = 0b00010000;

        const MADE_BROKEN_AND_ORDER_CHANGED = Self::MADE_BROKEN.bits() | Self::ORDER_CHANGED.bits();
        const CENTER_AND_MADE_BROKEN = Self::CENTER.bits() | Self::MADE_BROKEN.bits();
        const CENTER_AND_ORDER_CHANGED = Self::CENTER.bits() | Self::ORDER_CHANGED.bits();
        const CENTER_AND_MADE_BROKEN_AND_ORDER_CHANGED = Self::CENTER.bits() | Self::MADE_BROKEN.bits() | Self::ORDER_CHANGED.bits();
    }
}

impl Default for BondReactingCenter {
    fn default() -> Self {
        BondReactingCenter::UNMARKED
    }
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

/// SGroup type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupType {
    Superatom,     // SUP
    MultipleGroup, // MUL
    RepeatingUnit, // SRU
    Monomer,       // MON
    Mer,           // MER
    Copolymer,     // COP
    Crosslink,     // CRO
    Modification,  // MOD
    Graft,         // GRA
    Component,     // COM
    Mixture,       // MIX
    Formulation,   // FOR
    Data,          // DAT
    AnyPolymer,    // ANY
    Generic,       // GEN
}

/// SGroup subtype for polymers
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupSubtype {
    Alternating, // ALT
    Random,      // RAN
    Block,       // BLO
}

/// SGroup connectivity types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupConnectivity {
    HeadToHead,    // HH
    HeadToTail,    // HT
    EitherUnknown, // EU
}

/// SGroup multiplier term (variable or integer)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupMultiplierTerm {
    Variable(char),
    Integer(u32),
}

/// SGroup multiplier arithmetic operator
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupMultiplierOp {
    Add,
    Sub,
    Mul,
    Div,
}

/// SGroup multiplier for repeating unit properties
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupMultiplier {
    Single(SGroupMultiplierTerm),
    Expression {
        left: SGroupMultiplierTerm,
        op: SGroupMultiplierOp,
        right: SGroupMultiplierTerm,
    },
}

/// SGroup bracket coordinates
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SGroupBracketCoords {
    pub bracket1: (f64, f64),
    pub bracket2: (f64, f64),
}

/// SGroup connecting bond
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SGroupConnectingBond {
    pub bond_index: usize,
    pub bond_vector: (f64, f64),
}

/// SGroup bracket style
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupBracketStyle {
    #[default]
    Default, // 0 = default brackets
    Curved, // 1 = curved (parenthetic) brackets
}

/// SGroup data type
#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum SGroupDataType {
    Formatted,
    Numeric,
    #[default]
    Text,
}

/// SGroup data
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SGroupData {
    pub field_type: SGroupDataType,
    pub field_units: Option<String>,
    pub query_identifier: Option<String>,
    pub data_query_operator: Option<String>,
    pub data_content: Option<Vec<String>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupDataDisplayType {
    Attached,
    Detached,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupDataDisplayPlacement {
    Absolute,
    Relative,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupDataDisplayUnits {
    None,
    DisplayUnits,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupDataDisplayChars {
    All,
    Number(u32),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SGroupDataDisplay {
    pub coords: (f64, f64),
    pub display_type: SGroupDataDisplayType,
    pub display_placement: SGroupDataDisplayPlacement,
    pub display_units: SGroupDataDisplayUnits,
    pub display_chars: SGroupDataDisplayChars,
}

/// SGroup (Substance group)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SGroup {
    pub label: Option<u32>,
    pub subscript: Option<String>,
    pub group_type: SGroupType,
    pub group_subtype: Option<SGroupSubtype>,
    pub connectivity: Option<SGroupConnectivity>,
    pub expansion: bool,
    pub multiplier: Option<SGroupMultiplier>,
    pub atom_indices: Vec<usize>,
    pub bond_indices: Vec<usize>,
    pub parent_atom_indices: Option<Vec<usize>>,
    pub correspondence: Option<Vec<usize>>,
    pub connecting_bond: Option<SGroupConnectingBond>,
    pub bracket_coords: Option<SGroupBracketCoords>,
    pub hierarchy_parent: Option<usize>,
    pub component_number: Option<u32>,
    pub bracket_style: Option<SGroupBracketStyle>,
    pub data: BTreeMap<String, SGroupData>,
    pub display: Option<SGroupDataDisplay>,
}

impl SGroup {
    pub fn new(group_type: SGroupType) -> Self {
        Self {
            label: None,
            subscript: None,
            group_type,
            group_subtype: None,
            connectivity: None,
            expansion: false,
            multiplier: None,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
            parent_atom_indices: None,
            correspondence: None,
            connecting_bond: None,
            bracket_coords: None,
            hierarchy_parent: None,
            component_number: None,
            bracket_style: None,
            data: BTreeMap::new(),
            display: None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum RGroupOccurrence {
    Exactly(u8),
    Range(u8, u8),
    GreaterThan(u8),
    FewerThan(u8),
}

impl Default for RGroupOccurrence {
    fn default() -> Self {
        RGroupOccurrence::GreaterThan(0)
    }
}

impl RGroupOccurrence {
    pub fn contains(&self, count: u8) -> bool {
        match self {
            RGroupOccurrence::Exactly(n) => *n == count,
            RGroupOccurrence::Range(n, m) => count >= *n && count <= *m,
            RGroupOccurrence::GreaterThan(n) => count > *n,
            RGroupOccurrence::FewerThan(n) => count < *n,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RGroup {
    pub label: Option<u32>,
    pub dependent_label: Option<u32>,
    pub rgroup_or_h: bool,
    pub occurrence: Vec<RGroupOccurrence>,
}

impl Default for RGroup {
    fn default() -> Self {
        Self::new(None)
    }
}

impl RGroup {
    pub fn new(label: Option<u32>) -> Self {
        Self {
            label,
            dependent_label: None,
            rgroup_or_h: false,
            occurrence: vec![RGroupOccurrence::GreaterThan(0)],
        }
    }

    pub fn from_symbol_bytes(input: &[u8]) -> Option<Self> {
        debug_assert!(input.len() <= 3, "R-group symbol must be 1-3 characters");

        if input.is_empty() || input[0] != b'R' {
            None
        } else if input.len() == 1 || input.len() == 2 && input[1] == b'#' {
            Some(Self::new(None))
        } else {
            let num_str = &input[1..];
            if num_str.len() == 1 {
                if num_str[0] < b'0' || num_str[0] > b'9' {
                    None
                } else {
                    let label = (num_str[0] - b'0') as u32;
                    if label == 0 {
                        Some(Self::new(None))
                    } else {
                        Some(Self::new(Some(label)))
                    }
                }
            } else if num_str.len() == 2 {
                if num_str[0] < b'0' || num_str[0] > b'9' || num_str[1] < b'0' || num_str[1] > b'9'
                {
                    None
                } else {
                    let label = ((num_str[0] - b'0') * 10 + (num_str[1] - b'0')) as u32;
                    if label == 0 {
                        Some(Self::new(None))
                    } else {
                        Some(Self::new(Some(label)))
                    }
                }
            } else {
                None
            }
        }
    }

    pub fn from_symbol_str(input: &str) -> Option<Self> {
        Self::from_symbol_bytes(input.as_bytes())
    }
}

impl Display for RGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.label.is_some() {
            write!(f, "R{}", self.label.unwrap_or(0))
        } else {
            write!(f, "R#")
        }
    }
}

/// Convert element symbol to [u8; 2] key for alphabetical sorting
fn element_symbol_key(element: umol_data::Element) -> [u8; 2] {
    let symbol = element.symbol();
    let bytes = symbol.as_bytes();
    [
        bytes.first().copied().unwrap_or(0),
        bytes.get(1).copied().unwrap_or(0),
    ]
}

/// Format sum formula according to Hill notation
fn format_sum_formula(
    c_count: usize,
    h_count: usize,
    atom_counts: std::collections::BTreeMap<[u8; 2], (umol_data::Element, usize)>,
    charge: i32,
) -> String {
    let mut sum_formula = String::new();

    // Carbon first
    if c_count > 1 {
        sum_formula.push_str(&format!("C{}", c_count));
    } else if c_count == 1 {
        sum_formula.push('C');
    }

    // Hydrogen second
    if h_count > 1 {
        sum_formula.push_str(&format!("H{}", h_count));
    } else if h_count == 1 {
        sum_formula.push('H');
    }

    // Other elements alphabetically by symbol
    for (_, (element, count)) in atom_counts {
        if count > 1 {
            sum_formula.push_str(&format!("{}{}", element, count));
        } else {
            sum_formula.push_str(&element.to_string());
        }
    }

    // Charge at the end
    if charge != 0 {
        sum_formula.push_str(&format!("{:+}", charge));
    }

    sum_formula
}

/// Atom data supplied by parsers before conversion into IR atoms.
pub struct AtomData {
    pub element: Element,
    pub isotope: Option<u32>,
    pub charge: Option<i8>,
    pub hydrogen_count: Option<u8>,
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
                radical: None,
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
                radical: None,
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
            order: b.order,
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[test]
    fn test_molecule_default() {
        let mol = Molecule::default();
        assert!(mol.atoms.is_empty());
        assert!(mol.bonds.is_empty());
        assert!(mol.rings.is_empty());
        assert_eq!(mol.source_format, SourceFormat::UNKNOWN);
    }

    #[test]
    fn test_extended_molecule_from_molecule() {
        let mol = Molecule {
            atoms: vec![Atom::from_aliphatic_atom(Element::C)],
            bonds: vec![],
            rings: vec![],
            source_format: SourceFormat::SMILES,
        };
        let ext = ExtendedMolecule::from_molecule(mol.clone());
        assert_eq!(ext.atoms.len(), 1);
        assert_eq!(ext.source_format, SourceFormat::SMILES);
        assert!(ext.sgroups.is_empty());
        assert!(ext.rgroups.is_empty());
    }

    #[test]
    fn test_extended_molecule_to_basic() {
        let mut ext = ExtendedMolecule::default();
        ext.atoms.push(ExtendedAtom::new(AtomSymbol::Element(Element::O)));
        ext.source_format = SourceFormat::MOL;
        ext.sgroups.insert(0, SGroup::new(SGroupType::Superatom));

        let mol = ext.to_basic();
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(mol.source_format, SourceFormat::MOL);
    }

    #[test]
    fn test_sgroup_serialize() {
        let sgroup = SGroup::new(SGroupType::Superatom);
        let yaml = serde_yaml::to_string(&sgroup).expect("Failed to serialize SGroup to YAML");
        let deserialized: SGroup =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize SGroup from YAML");
        assert_eq!(sgroup, deserialized);
    }

    #[rstest]
    #[case(b"R", RGroup::new(None))]
    #[case(b"R#", RGroup::new(None))]
    #[case(b"R0", RGroup::new(None))]
    #[case(b"R1", RGroup::new(Some(1)))]
    #[case(b"R12", RGroup::new(Some(12)))]
    fn test_rgroup_from_symbol_bytes(#[case] input: &[u8], #[case] expected: RGroup) {
        let symbol = RGroup::from_symbol_bytes(input);
        assert!(symbol.is_some());
        assert_eq!(symbol.unwrap(), expected);
    }

    #[test]
    fn test_rgroup_serialize() {
        let rgroup = RGroup::new(Some(1));
        let yaml = serde_yaml::to_string(&rgroup).expect("Failed to serialize RGroup to YAML");
        let deserialized: RGroup =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize RGroup from YAML");
        assert_eq!(rgroup, deserialized);
    }

    #[test]
    fn test_bond_order_query() {
        assert!(!BondOrder::Single.is_query());
        assert!(!BondOrder::Double.is_query());
        assert!(BondOrder::SingleOrDouble.is_query());
        assert!(BondOrder::Any.is_query());
    }

    #[test]
    fn test_bond_order_extended() {
        assert!(!BondOrder::Single.is_extended());
        assert!(BondOrder::Zero.is_extended());
        assert!(BondOrder::Quadruple.is_extended());
        assert!(BondOrder::Quintuple.is_extended());
        assert!(BondOrder::Sextuple.is_extended());
    }

    #[test]
    fn test_radical_conversion() {
        assert_eq!(AtomRadical::from_unpaired_e(0), None);
        assert_eq!(AtomRadical::from_unpaired_e(1), Some(AtomRadical::Doublet));
        assert_eq!(AtomRadical::from_unpaired_e(2), Some(AtomRadical::Triplet));
        assert_eq!(AtomRadical::from_unpaired_e(5), Some(AtomRadical::Other(5)));

        assert_eq!(AtomRadical::Singlet.to_unpaired_e(), 0);
        assert_eq!(AtomRadical::Doublet.to_unpaired_e(), 1);
        assert_eq!(AtomRadical::Triplet.to_unpaired_e(), 2);
        assert_eq!(AtomRadical::Other(5).to_unpaired_e(), 5);
    }
}
