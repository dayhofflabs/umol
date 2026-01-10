//! Bond IR for Table IR.

use std::collections::HashMap;

use super::error::ConversionError;
use crate::span::Span;

/// Ordered pair of atom indices, `first <= second`
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AtomPair {
    first: u32,
    second: u32,
}

impl AtomPair {
    /// Create a new AtomPair, normalizing to ensure first <= second.
    pub fn new(a: u32, b: u32) -> Self {
        if a <= b {
            Self {
                first: a,
                second: b,
            }
        } else {
            Self {
                first: b,
                second: a,
            }
        }
    }

    /// Get the first (smaller) atom index.
    pub fn first(&self) -> u32 {
        self.first
    }

    /// Get the second (larger) atom index.
    pub fn second(&self) -> u32 {
        self.second
    }

    /// Get both atom indices as a tuple (first, second).
    pub fn as_tuple(&self) -> (u32, u32) {
        (self.first, self.second)
    }
}

// Basic Bond IR
#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    pub atoms: AtomPair,
    pub order: BondOrder,
    pub ring: Option<u32>,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDirection>,
    pub span: Option<Span>,
}

impl Bond {
    pub fn new(start_atom: u32, end_atom: u32, order: BondOrder) -> Self {
        Self {
            atoms: AtomPair::new(start_atom, end_atom),
            order,
            ring: None,
            stereo: None,
            direction: None,
            span: None,
        }
    }

    /// Create a bond with just the order (atom indices set to 0, for later update)
    pub fn with_order(order: BondOrder) -> Self {
        Self {
            atoms: AtomPair::new(0, 0),
            order,
            ring: None,
            stereo: None,
            direction: None,
            span: None,
        }
    }

    /// Get the start (first/smaller) atom index.
    pub fn start_atom(&self) -> u32 {
        self.atoms.first()
    }

    /// Get the end (second/larger) atom index.
    pub fn end_atom(&self) -> u32 {
        self.atoms.second()
    }

    /// Set the atom indices, normalizing automatically.
    pub fn set_atoms(&mut self, a: u32, b: u32) {
        self.atoms = AtomPair::new(a, b);
    }
}

/// Bond order
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Sextuple,
    Aromatic,
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
    Any,
}

impl BondOrder {
    pub fn from_value(value: u8) -> Option<Self> {
        match value {
            0 => Some(BondOrder::Zero),
            1 => Some(BondOrder::Single),
            2 => Some(BondOrder::Double),
            3 => Some(BondOrder::Triple),
            4 => Some(BondOrder::Quadruple),
            5 => Some(BondOrder::Quintuple),
            6 => Some(BondOrder::Sextuple),
            _ => None,
        }
    }

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
/// The wedge (pointed) end of the stereo bond is at the first atom
/// In SMILES: Up=/, Down=\
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondDirection {
    Up,        // MOL: Wedge (code 1), SMILES: /
    Down,      // MOL: Dash (code 6), SMILES: \
    Either,    // MOL code 4 (Either)
}

/// Double-bond stereochemistry (E/Z) annotation in IR
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondStereo {
    Cis,
    Trans,
    Either,
}

/// Extended bond IR
/// Temporary container for bond features of generalized molecules.
/// Includes extended bond fields from CTFile and SMILES parsers.
/// TODO: Split into multiple semantically defined structures.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedBond {
    pub atoms: AtomPair,
    pub order: BondOrder,
    pub ring: Option<u32>,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDirection>,
    pub topology: Option<BondTopology>,
    pub reacting_center: Option<BondReactingCenter>,
    pub properties: HashMap<String, String>,
    pub span: Option<Span>,
}
// End of TODO

impl ExtendedBond {
    pub fn new(start_atom: u32, end_atom: u32, order: BondOrder) -> Self {
        Self {
            atoms: AtomPair::new(start_atom, end_atom),
            order,
            ring: None,
            stereo: None,
            direction: None,
            topology: None,
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        }
    }

    /// Create a bond with just the order (atom indices set to 0, for later update)
    pub fn with_order(order: BondOrder) -> Self {
        Self {
            atoms: AtomPair::new(0, 0),
            order,
            ring: None,
            stereo: None,
            direction: None,
            topology: None,
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        }
    }

    /// Get the start (first/smaller) atom index.
    pub fn start_atom(&self) -> u32 {
        self.atoms.first()
    }

    /// Get the end (second/larger) atom index.
    pub fn end_atom(&self) -> u32 {
        self.atoms.second()
    }

    /// Set the atom indices, normalizing automatically.
    pub fn set_atoms(&mut self, a: u32, b: u32) {
        self.atoms = AtomPair::new(a, b);
    }

    /// Check if this bond has extended features that would be lost in conversion to basic Bond.
    pub fn has_extended_features(&self) -> bool {
        self.order.is_query()
            || self.order.is_extended()
            || self.topology.map_or(false, |t| !t.is_default())
            || self.reacting_center.map_or(false, |r| !r.is_default())
            || !self.properties.is_empty()
    }
}

impl From<Bond> for ExtendedBond {
    fn from(bond: Bond) -> Self {
        Self {
            atoms: bond.atoms,
            order: bond.order,
            ring: bond.ring,
            stereo: bond.stereo,
            direction: bond.direction,
            topology: None,
            reacting_center: None,
            properties: HashMap::new(),
            span: bond.span,
        }
    }
}

impl TryFrom<ExtendedBond> for Bond {
    type Error = ConversionError;

    fn try_from(extended: ExtendedBond) -> Result<Self, Self::Error> {
        // Check for extended features
        if extended.has_extended_features() {
            return Err(ConversionError::HasExtendedFeatures);
        }

        Ok(Self {
            atoms: extended.atoms,
            order: extended.order,
            ring: extended.ring,
            stereo: extended.stereo,
            direction: extended.direction,
            span: extended.span,
        })
    }
}

/// Bond topology (chain, ring, either) query
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondTopology {
    Chain,  // MOL code 2
    Ring,   // MOL code 1
    Either, // MOL code 0 (default/unspecified)
}

impl BondTopology {
    /// Returns true if this is the default (Either) topology
    pub fn is_default(&self) -> bool {
        matches!(self, BondTopology::Either)
    }
}

bitflags::bitflags! {
    /// Bond reacting center (from CTAB reactions) - bitflags
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct BondReactingCenter: u8 {
        const UNMARKED         = 0;
        const CENTER           = 1;
        const NOT_CENTER       = 1 << 1;
        const NO_CHANGE        = 1 << 2;
        const MADE_BROKEN      = 1 << 3;
        const ORDER_CHANGED    = 1 << 4;

        const MADE_BROKEN_AND_ORDER_CHANGED = Self::MADE_BROKEN.bits() | Self::ORDER_CHANGED.bits();
        const CENTER_AND_MADE_BROKEN = Self::CENTER.bits() | Self::MADE_BROKEN.bits();
        const CENTER_AND_ORDER_CHANGED = Self::CENTER.bits() | Self::ORDER_CHANGED.bits();
        const CENTER_AND_MADE_BROKEN_AND_ORDER_CHANGED = Self::CENTER.bits() | Self::MADE_BROKEN.bits() | Self::ORDER_CHANGED.bits();
    }
}

impl BondReactingCenter {
    /// Returns true if this is the default (UNMARKED) reacting center
    pub fn is_default(&self) -> bool {
        self.is_empty() || *self == BondReactingCenter::UNMARKED
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_pair_ordering() {
        let p1 = AtomPair::new(0, 1);
        let p2 = AtomPair::new(1, 0);
        assert_eq!(p1, p2);
        assert_eq!(p1.first(), 0);
        assert_eq!(p1.second(), 1);
    }

    #[test]
    fn test_atom_pair_ord() {
        let p1 = AtomPair::new(0, 1);
        let p2 = AtomPair::new(0, 2);
        let p3 = AtomPair::new(1, 2);
        assert!(p1 < p2);
        assert!(p2 < p3);
    }

    #[test]
    fn test_bond_new() {
        let bond = Bond::new(0, 1, BondOrder::Single);
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 1);
        assert_eq!(bond.order, BondOrder::Single);
    }

    #[test]
    fn test_bond_new_normalizes() {
        let bond = Bond::new(5, 2, BondOrder::Double);
        assert_eq!(bond.start_atom(), 2);
        assert_eq!(bond.end_atom(), 5);
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn test_bond_set_atoms() {
        let mut bond = Bond::with_order(BondOrder::Single);
        bond.set_atoms(7, 3);
        assert_eq!(bond.start_atom(), 3);
        assert_eq!(bond.end_atom(), 7);
    }

    #[test]
    fn test_bond_with_order() {
        let bond = Bond::with_order(BondOrder::Double);
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 0);
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn test_extended_bond_new() {
        let bond = ExtendedBond::new(0, 1, BondOrder::Triple);
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 1);
        assert_eq!(bond.order, BondOrder::Triple);
        assert!(!bond.has_extended_features());
    }

    #[test]
    fn test_extended_bond_new_normalizes() {
        let bond = ExtendedBond::new(10, 4, BondOrder::Aromatic);
        assert_eq!(bond.start_atom(), 4);
        assert_eq!(bond.end_atom(), 10);
        assert_eq!(bond.order, BondOrder::Aromatic);
    }

    #[test]
    fn test_extended_bond_set_atoms() {
        let mut bond = ExtendedBond::with_order(BondOrder::Double);
        bond.set_atoms(9, 1);
        assert_eq!(bond.start_atom(), 1);
        assert_eq!(bond.end_atom(), 9);
    }

    #[test]
    fn test_extended_bond_with_order() {
        let bond = ExtendedBond::with_order(BondOrder::Aromatic);
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 0);
        assert_eq!(bond.order, BondOrder::Aromatic);
        assert!(!bond.has_extended_features());
    }

    #[test]
    fn test_from_bond_to_extended_bond() {
        let bond = Bond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Single,
            ring: Some(5),
            stereo: Some(BondStereo::Cis),
            direction: Some(BondDirection::Up),
            span: None,
        };

        let extended: ExtendedBond = bond.into();
        assert_eq!(extended.start_atom(), 0);
        assert_eq!(extended.end_atom(), 1);
        assert_eq!(extended.order, BondOrder::Single);
        assert_eq!(extended.ring, Some(5));
        assert_eq!(extended.stereo, Some(BondStereo::Cis));
        assert_eq!(extended.direction, Some(BondDirection::Up));
        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_try_from_extended_bond_to_bond() {
        let extended = ExtendedBond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Double,
            ring: Some(6),
            stereo: Some(BondStereo::Trans),
            direction: Some(BondDirection::Down),
            topology: None,
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        };

        let bond: Bond = extended.try_into().unwrap();
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 1);
        assert_eq!(bond.order, BondOrder::Double);
        assert_eq!(bond.ring, Some(6));
        assert_eq!(bond.stereo, Some(BondStereo::Trans));
        assert_eq!(bond.direction, Some(BondDirection::Down));
    }

    #[test]
    fn test_try_from_extended_bond_to_bond_invalid() {
        let extended = ExtendedBond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Single,
            ring: None,
            stereo: None,
            direction: None,
            topology: Some(BondTopology::Ring),
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        };

        let result: Result<Bond, _> = extended.try_into();
        assert!(result.is_err());
    }

    #[test]
    fn test_has_extended_features_basic() {
        let extended = ExtendedBond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Triple,
            ring: Some(5),
            stereo: Some(BondStereo::Cis),
            direction: Some(BondDirection::Up),
            topology: None,
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        };

        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_has_extended_features_extended() {
        let extended = ExtendedBond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Single,
            ring: None,
            stereo: None,
            direction: None,
            topology: Some(BondTopology::Chain),
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        };

        assert!(extended.has_extended_features());
    }

    #[test]
    fn test_has_extended_features_query() {
        let extended = ExtendedBond::with_order(BondOrder::Any);
        assert!(extended.has_extended_features());
    }

    #[test]
    fn test_has_extended_features_extended_order() {
        let extended = ExtendedBond::with_order(BondOrder::Zero);
        assert!(extended.has_extended_features());
    }

    #[test]
    fn test_roundtrip_bond_to_extended_to_bond() {
        let bond = Bond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Aromatic,
            ring: Some(6),
            stereo: Some(BondStereo::Either),
            direction: Some(BondDirection::Either),
            span: None,
        };

        let extended: ExtendedBond = bond.clone().into();
        let bond2: Bond = extended.try_into().unwrap();

        assert_eq!(bond.start_atom(), bond2.start_atom());
        assert_eq!(bond.end_atom(), bond2.end_atom());
        assert_eq!(bond.order, bond2.order);
        assert_eq!(bond.ring, bond2.ring);
        assert_eq!(bond.stereo, bond2.stereo);
        assert_eq!(bond.direction, bond2.direction);
        assert_eq!(bond.span, bond2.span);
    }
}
