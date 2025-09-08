//! Intermediate representation for molecular structures

use std::collections::HashMap;

use umol::Result;
use umol_data::Element;

use crate::io::config::ParseConfig;
use crate::io::ctab::atom::{AtomRadical, AtomStereoParity};
use crate::io::ctab::bond::{BondReactingCenter, BondStereo, BondTopology};
use crate::io::ctab::sgroup::{SGroupConnectivity, SGroupDataType};

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
    pub fragments: Vec<Fragment>, // All structural units
    pub links: Vec<Link>,         // All connections between fragments

    // Property associations (chemical decomposition of data S-groups)
    pub property_annotations: Vec<PropertyAnnotation>,

    // Metadata
    pub header: Option<Header>,
    pub properties: HashMap<String, String>,
    pub source_format: SourceFormat,
    pub parsing_warnings: Vec<ParsingWarning>,
}

impl Molecule {
    /// Create a new empty molecule
    pub fn new(source_format: SourceFormat) -> Self {
        Self {
            source_format,
            ..Default::default()
        }
    }

    /// Check if the molecule contains variable substitution features
    pub fn is_variable(&self) -> bool {
        self.fragments.iter().any(|f| {
            matches!(
                f.fragment_type,
                FragmentType::CoreScaffold | FragmentType::VariablePlaceholder(_)
            )
        }) || self.atoms.iter().any(|a| a.rgroup_label.is_some())
    }

    /// Check if the molecule contains polymer features
    pub fn is_polymer(&self) -> bool {
        self.fragments
            .iter()
            .any(|f| matches!(f.fragment_type, FragmentType::RepeatingUnit { .. }))
            || self.links.iter().any(|l| {
                matches!(
                    l.link_type,
                    LinkType::PolymerChain { .. } | LinkType::CrossLink { .. }
                )
            })
    }

    /// Check if the molecule contains query features
    pub fn is_query(&self) -> bool {
        self.atoms.iter().any(|a| a.is_query())
            || self.bonds.iter().any(|b| b.is_query())
    }

    /// Check if the molecule has 3D coordinates
    pub fn is_3d(&self) -> bool {
        self.atoms.iter().any(|a| a.position.is_some())
    }

    /// Get all R-group labels used in the molecule
    pub fn rgroup_labels(&self) -> Vec<String> {
        self.atoms
            .iter()
            .filter_map(|a| a.rgroup_label.as_ref())
            .cloned()
            .collect()
    }
}

/// Atom IR
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    // Core atomic properties (common to all formats)
    pub element_or_query: ElementOrQuery,
    pub position: Option<Point3D>,
    pub charge: i8,
    pub isotope: Option<u32>,
    pub radical: Option<AtomRadical>,

    // Stereochemistry
    pub stereo_parity: Option<AtomStereoParity>,
    pub stereo_care: Option<bool>,
    pub inversion_retention: Option<InversionRetention>,

    // Query properties (MOL/SMARTS)
    pub hydrogen_count: Option<HydrogenConstraint>,
    pub valence: Option<ValenceConstraint>,
    pub ring_bond_count: Option<RingBondConstraint>,
    pub substitution_count: Option<SubstitutionCountConstraint>,
    pub unsaturated: Option<bool>,
    pub exact_change: Option<bool>,

    // Variable substitution
    pub rgroup_label: Option<String>, // "R1", "R2", etc.
    pub attachment_point: Option<AttachmentPointType>,
    pub attachment_order: Option<Vec<(usize, u8)>>,

    // Reaction mapping
    pub atom_map_num: Option<u32>,

    // SMILES-specific
    pub aromaticity_specified: Option<bool>,
    pub chirality_specified: Option<SmilesChirality>, // @ vs @@

    // Metadata
    pub source_format: SourceFormat,
    pub original_text: Option<String>,
}

impl Atom {
    /// Create a new empty atom
    pub fn new_element(element: Element, source_format: SourceFormat) -> Self {
        Self {
            element_or_query: ElementOrQuery::Element(element),
            position: None,
            charge: 0,
            isotope: None,
            radical: None,
            stereo_parity: None,
            stereo_care: None,
            inversion_retention: None,
            hydrogen_count: None,
            valence: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            exact_change: None,
            rgroup_label: None,
            attachment_point: None,
            attachment_order: None,
            atom_map_num: None,
            aromaticity_specified: None,
            chirality_specified: None,
            source_format,
            original_text: None,
        }
    }

    /// Check if this atom has any query features
    pub fn is_query(&self) -> bool {
        matches!(
            self.element_or_query,
            ElementOrQuery::QueryAtom(_) | ElementOrQuery::AtomList { .. }
        ) || self.hydrogen_count.is_some()
            || self.valence.is_some()
            || self.ring_bond_count.is_some()
            || self.substitution_count.is_some()
            || self.unsaturated.is_some()
    }

    /// Check if this atom is an R-group site
    pub fn is_rgroup(&self) -> bool {
        self.rgroup_label.is_some()
            || matches!(self.element_or_query, ElementOrQuery::VariationPoint(_))
    }
}

/// Bond IR
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    pub atom_indices: (usize, usize),

    // Core properties
    pub bond_type: BondTypeOrQuery,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDirection>,

    // Query properties
    pub topology: Option<BondTopology>,
    pub reacting_center: Option<BondReactingCenter>,

    // SMILES-specific (future)
    pub aromaticity_specified: Option<bool>, // Was this ':' bond?

    // Metadata
    pub source_format: SourceFormat,
    pub original_text: Option<String>,
}

impl Bond {
    /// Check if this bond has any query features
    pub fn is_query(&self) -> bool {
        matches!(
            self.bond_type,
            BondTypeOrQuery::SingleOrDouble
                | BondTypeOrQuery::SingleOrAromatic
                | BondTypeOrQuery::DoubleOrAromatic
                | BondTypeOrQuery::Any
                | BondTypeOrQuery::AromaticOrAliphatic
                | BondTypeOrQuery::RingBond
                | BondTypeOrQuery::NonRingBond
        ) || self.topology.is_some()
            || self.reacting_center.map_or(false, |rc| !rc.is_empty())
    }
}

/// Generalized structural unit (atoms + bonds + attachment sites)
#[derive(Debug, Clone, PartialEq)]
pub struct Fragment {
    pub atoms: Vec<usize>,                     // Atoms belonging to this fragment
    pub bonds: Vec<usize>,                     // Bonds belonging to this fragment
    pub attachment_sites: Vec<AttachmentSite>, // Structured connection points
    pub fragment_type: FragmentType,           // What kind of fragment this is
}

/// Structured connection point on a fragment
#[derive(Debug, Clone, PartialEq)]
pub struct AttachmentSite {
    pub atom_index: usize,                     // Which atom in fragment
    pub site_label: Option<String>,            // "R1", "head", "tail", etc.
    pub multiplicity: u8,                      // How many connections possible (1=mono, 2=bi)
    pub directionality: Option<SiteDirection>, // For oriented connections
    pub attachment_type: AttachmentPointType,  // First, Second, Both
}

/// Connection between fragments at their attachment sites
#[derive(Debug, Clone, PartialEq)]
pub struct Link {
    pub from_fragment: usize,                          // Source fragment ID
    pub from_site: usize,                              // Source attachment site ID
    pub to_fragment: usize,                            // Target fragment ID
    pub to_site: usize,                                // Target attachment site ID
    pub link_type: LinkType,                           // How they connect
    pub internal_structure: Option<Vec<usize>>,        // Bridge atoms if needed
    pub connectivity_rules: Option<ConnectivityRules>, // Occurrence, dependency logic
}

/// What kind of structural fragment this represents
#[derive(Debug, Clone, PartialEq)]
pub enum FragmentType {
    // Variable substitution chemistry
    CoreScaffold,                // Invariant template core
    VariablePlaceholder(String), // R1, R2, etc.

    // Macromolecular structure
    RepeatingUnit {
        repetition: RepetitionPattern,
        tacticity: Option<Tacticity>,
        subtype: Option<String>, // ALT, RAN, BLO from S-group subtype
    },

    // Structural recognition
    AbbreviatedStructure {
        abbreviation: String, // "Ph", "OMe", "tBu"
        subscript: Option<String>,
        expansion: bool,
    },

    // Chemical composition
    ComponentUnit {
        component_number: Option<u32>,
        is_mixture: bool,
    },
}

/// How fragments connect to each other
#[derive(Debug, Clone, PartialEq)]
pub enum LinkType {
    // Variable substitution
    SubstitutionBond, // R-group to core attachment

    // Polymer connections
    PolymerChain {
        connectivity: ConnectivityType,  // HeadToTail, etc.
        connection_type: ConnectionType, // LinearSequence, etc.
    },

    // Physical associations
    CrossLink {
        link_type: CrossLinkType, // Covalent, Ionic, etc.
        geometry: Option<CrossLinkGeometry>,
    },

    ComponentAssociation {
        connectivity: Option<SGroupConnectivity>, // HH/HT/EU
    },
}

/// Rules governing how links can be formed
#[derive(Debug, Clone, PartialEq)]
pub struct ConnectivityRules {
    pub occurrence: Vec<OccurrenceConstraint>, // How many times this can happen
    pub dependent_label: Option<String>,       // Depends on another fragment
    pub can_be_hydrogen: bool,                 // Whether H is valid
}

/// Directional orientation for attachment sites
#[derive(Debug, Clone, PartialEq)]
pub enum SiteDirection {
    Forward,       // 5' to 3', N to C terminus
    Reverse,       // 3' to 5', C to N terminus
    Bidirectional, // Can grow in either direction
}

// Additional enums needed for Fragment/Link types
#[derive(Debug, Clone, PartialEq)]
pub enum RepetitionPattern {
    Count(u32),       // Exact number
    Variable(String), // "n", "m", etc.
    Range(u32, u32),  // n to m repetitions
    Unlimited,        // Indefinite repetition
}

#[derive(Debug, Clone, PartialEq)]
pub enum Tacticity {
    Isotactic,    // All stereocenters same configuration
    Syndiotactic, // Alternating stereocenter configuration
    Atactic,      // Random stereocenter configuration
}

#[derive(Debug, Clone, PartialEq)]
pub enum ConnectivityType {
    HeadToTail,  // Normal polymer growth
    HeadToHead,  // Reverse connection
    TailToTail,  // Double reverse
    Random,      // No specific pattern
    Alternating, // A-B-A-B pattern
    Block,       // AAA-BBB pattern
}

#[derive(Debug, Clone, PartialEq)]
pub enum ConnectionType {
    LinearSequence, // End-to-end connection
    BranchPoint,    // Side chain attachment
    Crosslink,      // Inter-chain connection
}

#[derive(Debug, Clone, PartialEq)]
pub enum CrossLinkType {
    Covalent,    // Chemical bond
    Ionic,       // Electrostatic interaction
    Hydrogen,    // Hydrogen bonding
    VanDerWaals, // Weak interaction
}

#[derive(Debug, Clone, PartialEq)]
pub enum CrossLinkGeometry {
    Linear,              // Direct connection
    Bridged(Vec<usize>), // Connection through bridge atoms
}

/// Property data associated with molecular substructures (from Data S-groups)
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyAnnotation {
    pub target: AnnotationTarget,
    pub field_name: String,                  // Field identifier
    pub field_type: SGroupDataType,          // Formatted/Numeric/Text (direct from MOL)
    pub field_units: Option<String>,         // Units (direct from MOL)
    pub data_content: Vec<String>,           //  data strings (direct from MOL)
    pub query_identifier: Option<String>,    // Query ID (direct from MOL)
    pub data_query_operator: Option<String>, // Query operator (direct from MOL)
    pub display_info: Option<DisplayInfo>,   // Display coordinates, etc.
}

/// What the annotation refers to
#[derive(Debug, Clone, PartialEq)]
pub enum AnnotationTarget {
    Atom(usize),
    Bond(usize),
    Fragment(Vec<usize>), // Set of atoms/bonds
    Molecule,             // Whole molecule
}

/// Unified atom representation (elements + queries)
#[derive(Debug, Clone, PartialEq)]
pub enum ElementOrQuery {
    Element(Element),
    NamedIsotope(NamedIsotope),
    QueryAtom(QueryAtomType),
    AtomList {
        elements: Vec<Element>,
        exclude: bool,
    },
    VariationPoint(String), // R1, R2, X, Y, Ar, etc.
}

/// Extended query atom types (superset of MOL + SMARTS)
#[derive(Debug, Clone, PartialEq)]
pub enum QueryAtomType {
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
#[derive(Debug, Clone, PartialEq)]
pub enum BondTypeOrQuery {
    // Concrete bond types
    Single,
    Double,
    Triple,
    Aromatic,

    // Query bond types
    SingleOrDouble,   // MOL code 5
    SingleOrAromatic, // MOL code 6
    DoubleOrAromatic, // MOL code 7
    Any,              // MOL code 8
    Zero,             // Zero-order bond

    // SMILES/SMARTS extensions (future)
    AromaticOrAliphatic, // SMARTS ~
    RingBond,            // SMARTS @
    NonRingBond,         // SMARTS !@
}

/// Bond direction/wedging information
#[derive(Debug, Clone, PartialEq)]
pub enum BondDirection {
    Up,     // Wedge up from first atom
    Down,   // Dash down from first atom
    Either, // Unspecified direction
}

/// Hydrogen count constraints for query atoms
#[derive(Debug, Clone, PartialEq)]
pub enum HydrogenConstraint {
    Exact(u8),     // Exactly n hydrogens
    AtLeast(u8),   // At least n hydrogens
    AtMost(u8),    // At most n hydrogens
    Range(u8, u8), // Between n and m hydrogens
    None,          // No implicit hydrogens
}

/// Valence constraints for query atoms
#[derive(Debug, Clone, PartialEq)]
pub enum ValenceConstraint {
    Exact(u8),     // Exactly this valence
    List(Vec<u8>), // One of these valences
    Range(u8, u8), // Valence in this range
}

/// Ring bond count constraints
#[derive(Debug, Clone, PartialEq)]
pub enum RingBondConstraint {
    AsDrawn,     // As depicted in structure
    Exact(u8),   // Exactly n ring bonds
    AtLeast(u8), // At least n ring bonds
    AtMost(u8),  // At most n ring bonds
    None,        // No ring bonds
}

/// Substitution count constraints  
#[derive(Debug, Clone, PartialEq)]
pub enum SubstitutionCountConstraint {
    AsDrawn,     // As depicted in structure
    Exact(u8),   // Exactly n substituents
    AtLeast(u8), // At least n substituents
    AtMost(u8),  // At most n substituents
    None,        // No substituents
}

/// Bond type constraints for attachments
#[derive(Debug, Clone, PartialEq)]
pub enum BondTypeConstraint {
    Exact(BondTypeOrQuery),
    OneOf(Vec<BondTypeOrQuery>),
    AnyExcept(Vec<BondTypeOrQuery>),
}

/// Occurrence constraints (from RGroupOccurrence)
#[derive(Debug, Clone, PartialEq)]
pub enum OccurrenceConstraint {
    Exactly(u8),
    Range(u8, u8),   // Inclusive
    GreaterThan(u8), // Default is > 0
    FewerThan(u8),
}

/// Header information (format-agnostic)
#[derive(Debug, Clone, PartialEq)]
pub struct Header {
    pub title: Option<String>,
    pub program_info: Option<String>,
    pub comment: Option<String>,
    pub timestamp: Option<String>,
    pub format_version: Option<String>,
}

/// Display information for annotations
#[derive(Debug, Clone, PartialEq)]
pub struct DisplayInfo {
    pub coordinates: Option<(f64, f64)>,
    pub font_size: Option<u8>,
    pub color: Option<String>,
    pub style: Option<DisplayStyle>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum DisplayStyle {
    Attached,
    Detached,
    Overlaid,
}

/// Parsing warnings and diagnostics
#[derive(Debug, Clone, PartialEq)]
pub struct ParsingWarning {
    pub severity: WarningSeverity,
    pub message: String,
    pub location: Option<SourceLocation>,
    pub suggestion: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum WarningSeverity {
    Info,
    Warning,
    Error,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SourceLocation {
    pub line: usize,
    pub column: usize,
    pub context: Option<String>,
}

/// SMILES chirality specifications
#[derive(Debug, Clone, PartialEq)]
pub enum SmilesChirality {
    Clockwise,        // @
    CounterClockwise, // @@
    Unspecified,      // No chirality specified
}

/// Inversion/retention for reaction centers
#[derive(Debug, Clone, PartialEq)]
pub enum InversionRetention {
    Inverted, // Configuration inverts
    Retained, // Configuration preserved
}

/// Attachment point types for R-groups
#[derive(Debug, Clone, PartialEq)]
pub enum AttachmentPointType {
    First,  // Primary attachment
    Second, // Secondary attachment
    Both,   // Can use either
}

/// Named isotopes (D, T, etc.)
#[derive(Debug, Clone, PartialEq)]
pub struct NamedIsotope {
    pub element: Element,
    pub mass_number: u32,
    pub symbol: String, // "D", "T"
}

/// 3D coordinate type
pub type Point3D = (f64, f64, f64);

/// Simple condition for R-group dependencies
#[derive(Debug, Clone, PartialEq)]
pub enum SubstituentCondition {
    IsElement(Element),
    IsQuery(QueryAtomType),
    MatchesPattern(String), // SMARTS pattern
}

/// Trait for types that can be constructed from parsed molecular data
pub trait ParseTarget: Sized {
    fn allows_query_features() -> bool;
    fn allows_variable_substitution() -> bool;
    fn allows_polymer_features() -> bool;

    fn from_parsed_data(parsed: Molecule, config: &ParseConfig) -> Result<Self>;
}
