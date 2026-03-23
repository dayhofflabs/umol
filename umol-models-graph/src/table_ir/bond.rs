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
    pub unpaired_electrons: Option<u8>,
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
            unpaired_electrons: None,
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
            unpaired_electrons: None,
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
            unpaired_electrons: None,
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
    pub unpaired_electrons: Option<u8>,
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
            unpaired_electrons: None,
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
            unpaired_electrons: None,
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
            unpaired_electrons: None,
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
            unpaired_electrons: bond.unpaired_electrons,
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
            unpaired_electrons: extended.unpaired_electrons,
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
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::normal(0, 1, AtomPair::new(0, 1))]
    #[case::swapped(1, 0, AtomPair::new(0, 1))]
    #[case::equal(1, 1, AtomPair::new(1, 1))]
    fn test_atom_pair_ordering(#[case] a: u32, #[case] b: u32, #[case] expected: AtomPair) {
        assert_eq!(AtomPair::new(a, b), expected);
    }

    #[rstest]
    #[case::first(AtomPair::new(0, 2), AtomPair::new(1, 2), true)]
    #[case::second(AtomPair::new(0, 1), AtomPair::new(0, 2), true)]
    #[case::swapped(AtomPair::new(1, 0), AtomPair::new(0, 2), true)]
    fn test_atom_pair_ord(#[case] a: AtomPair, #[case] b: AtomPair, #[case] expected: bool) {
        assert_eq!(a < b, expected);
    }

    #[rstest]
    #[case::first(AtomPair::new(0, 1), 0, Some(1))]
    #[case::second(AtomPair::new(0, 1), 1, Some(0))]
    #[case::none(AtomPair::new(0, 1), 2, None)]
    fn test_atom_pair_other(
        #[case] pair: AtomPair,
        #[case] index: u32,
        #[case] expected: Option<u32>,
    ) {
        assert_eq!(pair.other(index), expected);
    }

    #[rstest]
    #[case::shared(BondDonation::Shared, BondDonation::Shared)]
    #[case::donating(BondDonation::Donating, BondDonation::Accepting)]
    #[case::accepting(BondDonation::Accepting, BondDonation::Donating)]
    #[case(BondDonation::Accepting, BondDonation::Donating)]
    fn test_bond_donation_flip(#[case] donation: BondDonation, #[case] expected: BondDonation) {
        assert_eq!(donation.flip(), expected);
    }

    #[rstest]
    #[case::normal(0, 1, BondOrder::Single, Bond::new(0, 1, BondOrder::Single))]
    #[case::swapped(1, 0, BondOrder::Single, Bond::new(0, 1, BondOrder::Single))]
    #[case::equal(1, 1, BondOrder::Single, Bond::new(1, 1, BondOrder::Single))]
    #[case::double(5, 2, BondOrder::Double, Bond::new(2, 5, BondOrder::Double))]
    fn test_bond_new(
        #[case] a: u32,
        #[case] b: u32,
        #[case] order: BondOrder,
        #[case] expected: Bond,
    ) {
        let bond = Bond::new(a, b, order);
        assert_eq!(bond, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::normal_donating(0, 4, BondOrder::Single, BondDonation::Donating, Some(BondDonation::Donating))]
    #[case::swapped_donating(4, 0, BondOrder::Single, BondDonation::Donating, Some(BondDonation::Accepting))]
    #[case::swapped_accepting(5, 2, BondOrder::Single, BondDonation::Accepting, Some(BondDonation::Donating))]
    fn test_bond_new_dative(
        #[case] a: u32,
        #[case] b: u32,
        #[case] order: BondOrder,
        #[case] donation: BondDonation,
        #[case] expected: Option<BondDonation>,
    ) {
        let bond = Bond::new_dative(a, b, order, donation);
        assert_eq!(bond.donation, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::normal_donating(0, 4, BondOrder::Single, BondDonation::Donating, Some(BondDonation::Donating))]
    #[case::swapped_donating(4, 0, BondOrder::Single, BondDonation::Donating, Some(BondDonation::Accepting))]
    #[case::swapped_accepting(5, 2, BondOrder::Single, BondDonation::Accepting, Some(BondDonation::Donating))]
    fn test_bond_update_atoms(
        #[case] a: u32,
        #[case] b: u32,
        #[case] order: BondOrder,
        #[case] donation: BondDonation,
        #[case] expected: Option<BondDonation>,
    ) {
        let bond = Bond::new_dative(a, b, order, donation);
        let updated = bond.update_atoms(1, 3);
        assert_eq!(updated.donation, expected);
    }

    #[rstest]
    #[case::normal(0, 1, BondOrder::Any, ExtendedBond::new(0, 1, BondOrder::Any))]
    #[case::swapped(1, 0, BondOrder::Any, ExtendedBond::new(0, 1, BondOrder::Any))]
    #[case::equal(1, 1, BondOrder::Any, ExtendedBond::new(1, 1, BondOrder::Any))]
    fn test_extended_bond_new(
        #[case] a: u32,
        #[case] b: u32,
        #[case] order: BondOrder,
        #[case] expected: ExtendedBond,
    ) {
        let bond = ExtendedBond::new(a, b, order);
        assert_eq!(bond, expected);
    }

    #[rstest]
    #[case::normal_donating(
        0,
        4,
        BondOrder::Single,
        BondDonation::Donating,
        Some(BondDonation::Donating)
    )]
    #[case::swapped_donating(
        4,
        0,
        BondOrder::Single,
        BondDonation::Donating,
        Some(BondDonation::Accepting)
    )]
    #[case::swapped_accepting(
        5,
        2,
        BondOrder::Single,
        BondDonation::Accepting,
        Some(BondDonation::Donating)
    )]
    fn test_extended_bond_new_dative(
        #[case] a: u32,
        #[case] b: u32,
        #[case] order: BondOrder,
        #[case] donation: BondDonation,
        #[case] expected: Option<BondDonation>,
    ) {
        let bond = ExtendedBond::new_dative(a, b, order, donation);
        assert_eq!(bond.donation, expected);
    }

    #[rstest]
    #[case::normal_donating(
        0,
        4,
        BondOrder::Single,
        BondDonation::Donating,
        Some(BondDonation::Donating)
    )]
    #[case::swapped_donating(
        4,
        0,
        BondOrder::Single,
        BondDonation::Donating,
        Some(BondDonation::Accepting)
    )]
    #[case::swapped_accepting(
        5,
        2,
        BondOrder::Single,
        BondDonation::Accepting,
        Some(BondDonation::Donating)
    )]
    fn test_extended_bond_update_atoms(
        #[case] a: u32,
        #[case] b: u32,
        #[case] order: BondOrder,
        #[case] donation: BondDonation,
        #[case] expected: Option<BondDonation>,
    ) {
        let bond = ExtendedBond::new_dative(a, b, order, donation);
        let updated = bond.update_atoms(1, 3);
        assert_eq!(updated.donation, expected);
    }
    #[rstest]
    #[case::normal(
        Bond::new(0, 1, BondOrder::Single),
        ExtendedBond::new(0, 1, BondOrder::Single)
    )]
    fn test_bond_into_extended_bond(#[case] bond: Bond, #[case] expected: ExtendedBond) {
        let extended: ExtendedBond = bond.into();
        assert_eq!(extended, expected);
    }

    #[rstest]
    #[case::normal(
        ExtendedBond::new(0, 1, BondOrder::Double),
        Bond::new(0, 1, BondOrder::Double)
    )]
    fn test_extended_bond_try_into_bond(#[case] extended: ExtendedBond, #[case] expected: Bond) {
        let bond: Bond = extended.try_into().unwrap();
        assert_eq!(bond, expected);
    }

    #[rstest]
    #[case::query(ExtendedBond::new(0, 1, BondOrder::Any))]
    #[case::topology(ExtendedBond { atoms: AtomPair::new(0, 1), order: BondOrder::Single, charge: None, unpaired_electrons: None,
                                    multiplicity: None, ring: None, stereo: None, wedge: None, donation: None, noncovalent: None,
                                    topology: Some(BondTopology::Ring), reacting_center: None, properties: HashMap::new(), span: None })]
    fn test_extended_bond_try_into_bond_error(#[case] extended: ExtendedBond) {
        let result: Result<Bond, _> = extended.try_into();
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            ConversionError::HasExtendedFeatures
        ));
    }

    #[rstest]
    #[case::normal(ExtendedBond::new(0, 1, BondOrder::Single), false)]
    #[case::query(ExtendedBond::new(0, 1, BondOrder::Any), true)]
    #[case::order_zero(ExtendedBond::new(0, 1, BondOrder::Zero), true)]
    #[case::topology(ExtendedBond { atoms: AtomPair::new(0, 1), order: BondOrder::Single, charge: None, unpaired_electrons: None,
                                    multiplicity: None, ring: None, stereo: None, wedge: None, donation: None, noncovalent: None,
                                    topology: Some(BondTopology::Ring), reacting_center: None, properties: HashMap::new(), span: None }, true)]
    fn test_extended_bond_has_extended_features(
        #[case] extended: ExtendedBond,
        #[case] expected: bool,
    ) {
        assert_eq!(extended.has_extended_features(), expected);
    }
}
