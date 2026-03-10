//! Bond types for Table IR.

use std::collections::HashMap;

use strum::{Display, EnumString};
use umol_data::SpinMultiplicity;

use super::error::ConversionError;
use crate::bond::BondNoncovalent;
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

    /// Check if this bond contains the given atom index.
    pub fn contains(&self, index: u32) -> bool {
        self.first == index || self.second == index
    }

    /// Get the other atom index.
    pub fn other(&self, index: u32) -> Option<u32> {
        if self.first == index {
            Some(self.second)
        } else if self.second == index {
            Some(self.first)
        } else {
            None
        }
    }
}

/// Electron pair donation for dative/coordinate bonds.
/// Direction is defined from the perspective of the first (smaller-indexed) atom.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondDonation {
    Shared,    // Normal covalent bond (both atoms contribute)
    Donating,  // First atom donates electron pair to second
    Accepting, // First atom accepts electron pair from second
}

impl BondDonation {
    /// Flip the donation
    pub fn flip(self) -> Self {
        match self {
            Self::Shared => Self::Shared,
            Self::Donating => Self::Accepting,
            Self::Accepting => Self::Donating,
        }
    }
}

// Basic Bond IR
#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    pub atoms: AtomPair,
    pub order: BondOrder,
    pub charge: Option<i8>,
    pub multiplicity: Option<SpinMultiplicity>,
    pub ring: Option<u32>,
    pub stereo: Option<BondStereo>,
    pub wedge: Option<BondWedge>,
    pub donation: Option<BondDonation>,
    pub noncovalent: Option<BondNoncovalent>,
    pub span: Option<Span>,
}

impl Bond {
    pub fn new(a: u32, b: u32, order: BondOrder) -> Self {
        Self {
            atoms: AtomPair::new(a, b),
            order,
            charge: None,
            multiplicity: None,
            ring: None,
            stereo: None,
            wedge: None,
            donation: None,
            noncovalent: None,
            span: None,
        }
    }

    /// Create a dative bond, adjusting donation for AtomPair normalization.
    /// The donation parameter describes the donation from `a` to `b` before normalization.
    pub fn new_dative(a: u32, b: u32, order: BondOrder, donation: BondDonation) -> Self {
        let swapped = a > b;
        Self {
            atoms: AtomPair::new(a, b),
            order,
            charge: None,
            multiplicity: None,
            ring: None,
            stereo: None,
            wedge: None,
            donation: Some(if swapped { donation.flip() } else { donation }),
            noncovalent: None,
            span: None,
        }
    }

    /// Create a non-covalent bond (hydrogen bond, halogen bond, etc.)
    pub fn new_noncovalent(a: u32, b: u32, noncovalent: BondNoncovalent) -> Self {
        Self {
            atoms: AtomPair::new(a, b),
            order: BondOrder::Zero,
            charge: None,
            multiplicity: None,
            ring: None,
            stereo: None,
            wedge: None,
            donation: None,
            noncovalent: Some(noncovalent),
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

    /// Return a bond with updated atom indices while preserving donation semantics.
    pub fn update_atoms(&self, a: u32, b: u32) -> Self {
        let mut updated = self.clone();
        updated.atoms = AtomPair::new(a, b);
        if a > b {
            updated.donation = updated.donation.map(|d| d.flip());
        }
        updated
    }
}

/// Bond order
#[derive(Clone, Copy, Debug, PartialEq, Display, EnumString)]
#[strum(serialize_all = "snake_case")]
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
    /// Get the value of the bond order if well-defined
    /// Returns None for query bonds and aromatic bonds.
    pub fn value(&self) -> Option<u8> {
        match self {
            BondOrder::Zero => Some(0),
            BondOrder::Single => Some(1),
            BondOrder::Double => Some(2),
            BondOrder::Triple => Some(3),
            BondOrder::Quadruple => Some(4),
            BondOrder::Quintuple => Some(5),
            BondOrder::Sextuple => Some(6),
            BondOrder::Aromatic => None,
            BondOrder::SingleOrDouble => None,
            BondOrder::SingleOrAromatic => None,
            BondOrder::DoubleOrAromatic => None,
            BondOrder::Any => None,
        }
    }

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

/// Stereo wedge indication for single bonds
/// In MOL files: Wedge (up, code 1), Dash (down, code 6)
/// The wedge (pointed) end of the stereo bond is at the first atom
/// In SMILES: / (up), \ (down)
/// In CXSMILES: w: (undefined), wU: (undefined, display up), wD: (undefined, display down)
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondWedge {
    Up,         // MOL: Wedge (code 1), SMILES: /
    Down,       // MOL: Dash (code 6), SMILES: \
    Either,     // MOL code 4, CXSMILES: w: (stereo undefined)
    EitherUp,   // CXSMILES: wU: (stereo undefined, display up)
    EitherDown, // CXSMILES: wD: (stereo undefined, display down)
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
    pub charge: Option<i8>,
    pub multiplicity: Option<SpinMultiplicity>,
    pub ring: Option<u32>,
    pub stereo: Option<BondStereo>,
    pub wedge: Option<BondWedge>,
    pub donation: Option<BondDonation>,
    pub noncovalent: Option<BondNoncovalent>,
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
            charge: None,
            multiplicity: None,
            ring: None,
            stereo: None,
            wedge: None,
            donation: None,
            noncovalent: None,
            topology: None,
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        }
    }

    /// Create a dative bond, adjusting donation for AtomPair normalization.
    /// The donation parameter describes the donation from `a` to `b` before normalization.
    pub fn new_dative(a: u32, b: u32, order: BondOrder, donation: BondDonation) -> Self {
        let swapped = a > b;
        Self {
            atoms: AtomPair::new(a, b),
            order,
            charge: None,
            multiplicity: None,
            ring: None,
            stereo: None,
            wedge: None,
            donation: Some(if swapped { donation.flip() } else { donation }),
            noncovalent: None,
            topology: None,
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        }
    }

    /// Create a non-covalent bond (hydrogen bond, halogen bond, etc.)
    pub fn new_noncovalent(a: u32, b: u32, noncovalent: BondNoncovalent) -> Self {
        Self {
            atoms: AtomPair::new(a, b),
            order: BondOrder::Zero,
            charge: None,
            multiplicity: None,
            ring: None,
            stereo: None,
            wedge: None,
            donation: None,
            noncovalent: Some(noncovalent),
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

    /// Return a bond with updated atom indices while preserving donation semantics.
    pub fn update_atoms(&self, a: u32, b: u32) -> Self {
        let mut updated = self.clone();
        updated.atoms = AtomPair::new(a, b);
        if a > b {
            updated.donation = updated.donation.map(|d| d.flip());
        }
        updated
    }

    /// Check if this bond has extended features that would be lost in conversion to basic Bond.
    pub fn has_extended_features(&self) -> bool {
        self.order.is_query()
            || self.order.is_extended()
            || self.topology.is_some_and(|t| !t.is_default())
            || self.reacting_center.is_some_and(|r| !r.is_default())
            || !self.properties.is_empty()
    }
}

impl From<Bond> for ExtendedBond {
    fn from(bond: Bond) -> Self {
        Self {
            atoms: bond.atoms,
            order: bond.order,
            charge: bond.charge,
            multiplicity: bond.multiplicity,
            ring: bond.ring,
            stereo: bond.stereo,
            wedge: bond.wedge,
            donation: bond.donation,
            noncovalent: bond.noncovalent,
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
        if extended.has_extended_features() {
            return Err(ConversionError::HasExtendedFeatures);
        }

        Ok(Self {
            atoms: extended.atoms,
            order: extended.order,
            charge: extended.charge,
            multiplicity: extended.multiplicity,
            ring: extended.ring,
            stereo: extended.stereo,
            wedge: extended.wedge,
            donation: extended.donation,
            noncovalent: extended.noncovalent,
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
    fn test_atom_pair_other() {
        let p1 = AtomPair::new(0, 1);
        assert_eq!(p1.other(0), Some(1));
        assert_eq!(p1.other(1), Some(0));
        assert_eq!(p1.other(2), None);
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
    fn test_bond_donation_flip() {
        assert_eq!(BondDonation::Shared.flip(), BondDonation::Shared);
        assert_eq!(BondDonation::Donating.flip(), BondDonation::Accepting);
        assert_eq!(BondDonation::Accepting.flip(), BondDonation::Donating);
    }

    #[test]
    fn test_bond_new_dative_no_swap() {
        // a < b, no swap needed
        let bond = Bond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating);
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 4);
        assert_eq!(bond.donation, Some(BondDonation::Donating));
    }

    #[test]
    fn test_bond_new_dative_with_swap() {
        // a > b, swap occurs, donation should flip
        let bond = Bond::new_dative(4, 0, BondOrder::Single, BondDonation::Donating);
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 4);
        // Original: 4 donates to 0. After swap: 0 accepts from 4.
        assert_eq!(bond.donation, Some(BondDonation::Accepting));
    }

    #[test]
    fn test_bond_new_dative_accepting_with_swap() {
        // a > b, swap occurs
        let bond = Bond::new_dative(5, 2, BondOrder::Single, BondDonation::Accepting);
        assert_eq!(bond.start_atom(), 2);
        assert_eq!(bond.end_atom(), 5);
        // Original: 5 accepts from 2. After swap: 2 donates to 5.
        assert_eq!(bond.donation, Some(BondDonation::Donating));
    }

    #[test]
    fn test_bond_update_atoms_no_swap() {
        let bond = Bond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating);
        let updated = bond.update_atoms(1, 3);
        assert_eq!(updated.start_atom(), 1);
        assert_eq!(updated.end_atom(), 3);
        assert_eq!(updated.donation, Some(BondDonation::Donating));
    }

    #[test]
    fn test_bond_update_atoms_with_swap() {
        let bond = Bond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating);
        let updated = bond.update_atoms(3, 1);
        assert_eq!(updated.start_atom(), 1);
        assert_eq!(updated.end_atom(), 3);
        assert_eq!(updated.donation, Some(BondDonation::Accepting));
    }

    #[test]
    fn test_extended_bond_new_dative_no_swap() {
        let bond = ExtendedBond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating);
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 4);
        assert_eq!(bond.donation, Some(BondDonation::Donating));
    }

    #[test]
    fn test_extended_bond_new_dative_with_swap() {
        let bond = ExtendedBond::new_dative(4, 0, BondOrder::Single, BondDonation::Accepting);
        assert_eq!(bond.start_atom(), 0);
        assert_eq!(bond.end_atom(), 4);
        // Original: 4 accepts from 0. After swap: 0 donates to 4.
        assert_eq!(bond.donation, Some(BondDonation::Donating));
    }

    #[test]
    fn test_extended_bond_update_atoms_no_swap() {
        let bond = ExtendedBond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating);
        let updated = bond.update_atoms(1, 3);
        assert_eq!(updated.start_atom(), 1);
        assert_eq!(updated.end_atom(), 3);
        assert_eq!(updated.donation, Some(BondDonation::Donating));
    }

    #[test]
    fn test_extended_bond_update_atoms_with_swap() {
        let bond = ExtendedBond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating);
        let updated = bond.update_atoms(3, 1);
        assert_eq!(updated.start_atom(), 1);
        assert_eq!(updated.end_atom(), 3);
        assert_eq!(updated.donation, Some(BondDonation::Accepting));
    }

    #[test]
    fn test_from_bond_to_extended_bond() {
        let bond = Bond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Single,
            charge: None,
            multiplicity: None,
            ring: Some(5),
            stereo: Some(BondStereo::Cis),
            wedge: Some(BondWedge::Up),
            donation: None,
            noncovalent: None,
            span: None,
        };

        let extended: ExtendedBond = bond.into();
        assert_eq!(extended.start_atom(), 0);
        assert_eq!(extended.end_atom(), 1);
        assert_eq!(extended.order, BondOrder::Single);
        assert_eq!(extended.ring, Some(5));
        assert_eq!(extended.stereo, Some(BondStereo::Cis));
        assert_eq!(extended.wedge, Some(BondWedge::Up));
        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_try_from_extended_bond_to_bond() {
        let extended = ExtendedBond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Double,
            charge: None,
            multiplicity: None,
            ring: Some(6),
            stereo: Some(BondStereo::Trans),
            wedge: Some(BondWedge::Down),
            donation: None,
            noncovalent: None,
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
        assert_eq!(bond.wedge, Some(BondWedge::Down));
    }

    #[test]
    fn test_try_from_extended_bond_to_bond_invalid() {
        let extended = ExtendedBond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Single,
            charge: None,
            multiplicity: None,
            ring: None,
            stereo: None,
            wedge: None,
            donation: None,
            noncovalent: None,
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
            charge: None,
            multiplicity: None,
            ring: Some(5),
            stereo: Some(BondStereo::Cis),
            wedge: Some(BondWedge::Up),
            donation: None,
            noncovalent: None,
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
            charge: None,
            multiplicity: None,
            ring: None,
            stereo: None,
            wedge: None,
            donation: None,
            noncovalent: None,
            topology: Some(BondTopology::Chain),
            reacting_center: None,
            properties: HashMap::new(),
            span: None,
        };

        assert!(extended.has_extended_features());
    }

    #[test]
    fn test_has_extended_features_query() {
        let extended = ExtendedBond::new(0, 1, BondOrder::Any);
        assert!(extended.has_extended_features());
    }

    #[test]
    fn test_has_extended_features_extended_order() {
        let extended = ExtendedBond::new(0, 1, BondOrder::Zero);
        assert!(extended.has_extended_features());
    }

    #[test]
    fn test_roundtrip_bond_to_extended_to_bond() {
        let bond = Bond {
            atoms: AtomPair::new(0, 1),
            order: BondOrder::Aromatic,
            charge: None,
            multiplicity: None,
            ring: Some(6),
            stereo: Some(BondStereo::Either),
            wedge: Some(BondWedge::Either),
            donation: None,
            noncovalent: None,
            span: None,
        };

        let extended: ExtendedBond = bond.clone().into();
        let bond2: Bond = extended.try_into().unwrap();

        assert_eq!(bond.start_atom(), bond2.start_atom());
        assert_eq!(bond.end_atom(), bond2.end_atom());
        assert_eq!(bond.order, bond2.order);
        assert_eq!(bond.ring, bond2.ring);
        assert_eq!(bond.stereo, bond2.stereo);
        assert_eq!(bond.wedge, bond2.wedge);
        assert_eq!(bond.span, bond2.span);
    }
}
