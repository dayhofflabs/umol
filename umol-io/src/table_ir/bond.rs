//! Bond types for Table IR.

use std::collections::HashMap;

use strum::{Display, EnumString};
use umol_chem::spin::SpinMultiplicity;

use super::error::ConversionError;
use super::span::Span;

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

/// Non-covalent interaction type for weak bonds (H-bonds, halogen bonds, etc.)
/// These bonds do not contribute to valence calculations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondNoncovalent {
    Hydrogen,
}

// Basic Bond IR
#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    pub atoms: AtomPair,
    pub order: BondOrder,
    pub donation: Option<BondDonation>,
    pub noncovalent: Option<BondNoncovalent>,
    pub charge: Option<i8>,
    pub unpaired_electrons: Option<u8>,
    pub multiplicity: Option<SpinMultiplicity>,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDirection>,
    pub wedge: Option<BondWedge>,
    pub ring: Option<u32>,
    pub span: Option<Span>,
}

impl Bond {
    pub fn new(a: u32, b: u32, order: BondOrder) -> Self {
        Self {
            atoms: AtomPair::new(a, b),
            order,
            donation: None,
            noncovalent: None,
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
            stereo: None,
            direction: None,
            wedge: None,
            ring: None,
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
            donation: Some(if swapped { donation.flip() } else { donation }),
            noncovalent: None,
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
            stereo: None,
            direction: None,
            wedge: None,
            ring: None,
            span: None,
        }
    }

    /// Create a non-covalent bond (hydrogen bond, halogen bond, etc.)
    pub fn new_noncovalent(a: u32, b: u32, noncovalent: BondNoncovalent) -> Self {
        Self {
            atoms: AtomPair::new(a, b),
            order: BondOrder::Zero,
            donation: None,
            noncovalent: Some(noncovalent),
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
            stereo: None,
            direction: None,
            wedge: None,
            ring: None,
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

    /// Return the bond with `a` and `b` as the new indices of its first and second atoms.
    /// Reversing their order re-reads the donation and the wedge from the other endpoint.
    pub fn update_atoms(&self, a: u32, b: u32) -> Self {
        let mut updated = self.clone();
        updated.atoms = AtomPair::new(a, b);
        if a > b {
            updated.donation = updated.donation.map(BondDonation::flip);
            updated.wedge = updated.wedge.map(BondWedge::flip);
        }
        updated
    }

    /// Narrow (pointed) endpoint of the wedge, if the bond carries one.
    pub fn narrow_endpoint(&self) -> Option<u32> {
        self.wedge.map(|wedge| match wedge.taper {
            BondTaper::Widening => self.atoms.first(),
            BondTaper::Narrowing => self.atoms.second(),
        })
    }

    /// Wide endpoint of the wedge, if the bond carries one.
    pub fn wide_endpoint(&self) -> Option<u32> {
        self.wedge.map(|wedge| match wedge.taper {
            BondTaper::Widening => self.atoms.second(),
            BondTaper::Narrowing => self.atoms.first(),
        })
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

/// Out-of-plane reading of a tetrahedral depiction wedge on a single bond.
/// In MOL files: Up (code 1), Down (code 6), Either (code 4).
/// In CXSMILES: w: (Either), wU: (EitherUp), wD: (EitherDown).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondOrientation {
    Up,         // MOL: Wedge (code 1)
    Down,       // MOL: Dash (code 6)
    Either,     // MOL code 4, CXSMILES: w: (stereo undefined)
    EitherUp,   // CXSMILES: wU: (stereo undefined, display up)
    EitherDown, // CXSMILES: wD: (stereo undefined, display down)
}

/// Width of a wedge from the stored pair's first endpoint toward its second.
/// `Widening` puts the narrow (pointed) end at `first()`, `Narrowing` at `second()`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BondTaper {
    Widening,
    Narrowing,
}

impl BondTaper {
    /// Reverse the taper, as when the same bond is read from its other end.
    pub fn flip(self) -> Self {
        match self {
            Self::Widening => Self::Narrowing,
            Self::Narrowing => Self::Widening,
        }
    }
}

/// Depiction wedge of a bond: its out-of-plane reading and which stored endpoint is the narrow
/// (pointed) end. A wedge describes the configuration at its narrow endpoint only.
/// In MOL files the bond line's first atom is the narrow end; CXSMILES wiggly-bond entries name
/// the narrow-end atom explicitly.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct BondWedge {
    pub orientation: BondOrientation,
    pub taper: BondTaper,
}

impl BondWedge {
    /// The wedge as read from the bond's other end: the taper reverses, the orientation is
    /// unchanged.
    pub fn flip(self) -> Self {
        Self {
            orientation: self.orientation,
            taper: self.taper.flip(),
        }
    }
}

/// Cis/trans directional bond, for a single bond adjacent to a double bond.
/// In SMILES/CXSMILES: `/` (Rising), `\` (Falling).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BondDirection {
    Rising,
    Falling,
}

impl BondDirection {
    /// Reverse the direction, as when the same bond is read from its other end.
    pub fn flip(self) -> Self {
        match self {
            Self::Rising => Self::Falling,
            Self::Falling => Self::Rising,
        }
    }
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
    pub topology: Option<BondTopology>,
    pub donation: Option<BondDonation>,
    pub noncovalent: Option<BondNoncovalent>,
    pub charge: Option<i8>,
    pub unpaired_electrons: Option<u8>,
    pub multiplicity: Option<SpinMultiplicity>,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDirection>,
    pub wedge: Option<BondWedge>,
    pub ring: Option<u32>,
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
            topology: None,
            donation: None,
            noncovalent: None,
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
            stereo: None,
            direction: None,
            wedge: None,
            ring: None,
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
            topology: None,
            donation: Some(if swapped { donation.flip() } else { donation }),
            noncovalent: None,
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
            stereo: None,
            direction: None,
            wedge: None,
            ring: None,
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
            topology: None,
            donation: None,
            noncovalent: Some(noncovalent),
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
            stereo: None,
            direction: None,
            wedge: None,
            ring: None,
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

    /// Return the bond with `a` and `b` as the new indices of its first and second atoms.
    /// Reversing their order re-reads the donation and the wedge from the other endpoint.
    pub fn update_atoms(&self, a: u32, b: u32) -> Self {
        let mut updated = self.clone();
        updated.atoms = AtomPair::new(a, b);
        if a > b {
            updated.donation = updated.donation.map(BondDonation::flip);
            updated.wedge = updated.wedge.map(BondWedge::flip);
        }
        updated
    }

    /// Narrow (pointed) endpoint of the wedge, if the bond carries one.
    pub fn narrow_endpoint(&self) -> Option<u32> {
        self.wedge.map(|wedge| match wedge.taper {
            BondTaper::Widening => self.atoms.first(),
            BondTaper::Narrowing => self.atoms.second(),
        })
    }

    /// Wide endpoint of the wedge, if the bond carries one.
    pub fn wide_endpoint(&self) -> Option<u32> {
        self.wedge.map(|wedge| match wedge.taper {
            BondTaper::Widening => self.atoms.second(),
            BondTaper::Narrowing => self.atoms.first(),
        })
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
            topology: None,
            donation: bond.donation,
            noncovalent: bond.noncovalent,
            charge: bond.charge,
            unpaired_electrons: bond.unpaired_electrons,
            multiplicity: bond.multiplicity,
            stereo: bond.stereo,
            direction: bond.direction,
            wedge: bond.wedge,
            ring: bond.ring,
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
            donation: extended.donation,
            noncovalent: extended.noncovalent,
            charge: extended.charge,
            unpaired_electrons: extended.unpaired_electrons,
            multiplicity: extended.multiplicity,
            stereo: extended.stereo,
            direction: extended.direction,
            wedge: extended.wedge,
            ring: extended.ring,
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
    #[case::renumbered_donating(Bond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating), 1, 3, Bond::new_dative(1, 3, BondOrder::Single, BondDonation::Donating))]
    #[case::reversed_donating(Bond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating), 3, 1, Bond::new_dative(1, 3, BondOrder::Single, BondDonation::Accepting))]
    #[case::reversed_accepting(Bond::new_dative(2, 5, BondOrder::Single, BondDonation::Accepting), 3, 1, Bond::new_dative(1, 3, BondOrder::Single, BondDonation::Donating))]
    #[case::renumbered_wedge(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..Bond::new(0, 4, BondOrder::Single) }, 1, 3, Bond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..Bond::new(1, 3, BondOrder::Single) })]
    #[case::reversed_wedge(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..Bond::new(0, 4, BondOrder::Single) }, 3, 1, Bond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Narrowing }), ..Bond::new(1, 3, BondOrder::Single) })]
    #[case::reversed_wedge_and_donation(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }), ..Bond::new_dative(0, 4, BondOrder::Single, BondDonation::Accepting) }, 3, 1, Bond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Widening }), ..Bond::new_dative(1, 3, BondOrder::Single, BondDonation::Donating) })]
    fn test_bond_update_atoms(
        #[case] bond: Bond,
        #[case] a: u32,
        #[case] b: u32,
        #[case] expected: Bond,
    ) {
        assert_eq!(bond.update_atoms(a, b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::plain(Bond::new(0, 4, BondOrder::Single))]
    #[case::donating(Bond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating))]
    #[case::narrowing_wedge(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }), ..Bond::new(0, 4, BondOrder::Single) })]
    fn test_bond_update_atoms_identity(#[case] bond: Bond) {
        assert_eq!(bond.update_atoms(bond.start_atom(), bond.end_atom()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::none(Bond::new(0, 4, BondOrder::Single), None)]
    #[case::widening(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..Bond::new(0, 4, BondOrder::Single) }, Some(0))]
    #[case::narrowing(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Narrowing }), ..Bond::new(0, 4, BondOrder::Single) }, Some(4))]
    fn test_bond_narrow_endpoint(#[case] bond: Bond, #[case] expected: Option<u32>) {
        assert_eq!(bond.narrow_endpoint(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::none(Bond::new(0, 4, BondOrder::Single), None)]
    #[case::widening(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..Bond::new(0, 4, BondOrder::Single) }, Some(4))]
    #[case::narrowing(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Narrowing }), ..Bond::new(0, 4, BondOrder::Single) }, Some(0))]
    fn test_bond_wide_endpoint(#[case] bond: Bond, #[case] expected: Option<u32>) {
        assert_eq!(bond.wide_endpoint(), expected);
    }

    #[rstest]
    #[case::widening(BondTaper::Widening, BondTaper::Narrowing)]
    #[case::narrowing(BondTaper::Narrowing, BondTaper::Widening)]
    fn test_bond_taper_flip(#[case] taper: BondTaper, #[case] expected: BondTaper) {
        assert_eq!(taper.flip(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::up_widening(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }, BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Narrowing })]
    #[case::down_narrowing(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }, BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Widening })]
    #[case::either_widening(BondWedge { orientation: BondOrientation::Either, taper: BondTaper::Widening }, BondWedge { orientation: BondOrientation::Either, taper: BondTaper::Narrowing })]
    fn test_bond_wedge_flip(#[case] wedge: BondWedge, #[case] expected: BondWedge) {
        assert_eq!(wedge.flip(), expected);
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

    #[rustfmt::skip]
    #[rstest]
    #[case::renumbered_donating(ExtendedBond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating), 1, 3, ExtendedBond::new_dative(1, 3, BondOrder::Single, BondDonation::Donating))]
    #[case::reversed_donating(ExtendedBond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating), 3, 1, ExtendedBond::new_dative(1, 3, BondOrder::Single, BondDonation::Accepting))]
    #[case::reversed_accepting(ExtendedBond::new_dative(2, 5, BondOrder::Single, BondDonation::Accepting), 3, 1, ExtendedBond::new_dative(1, 3, BondOrder::Single, BondDonation::Donating))]
    #[case::renumbered_wedge(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..ExtendedBond::new(0, 4, BondOrder::Single) }, 1, 3, ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..ExtendedBond::new(1, 3, BondOrder::Single) })]
    #[case::reversed_wedge(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..ExtendedBond::new(0, 4, BondOrder::Single) }, 3, 1, ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Narrowing }), ..ExtendedBond::new(1, 3, BondOrder::Single) })]
    #[case::reversed_wedge_and_donation(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }), ..ExtendedBond::new_dative(0, 4, BondOrder::Single, BondDonation::Accepting) }, 3, 1, ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Widening }), ..ExtendedBond::new_dative(1, 3, BondOrder::Single, BondDonation::Donating) })]
    fn test_extended_bond_update_atoms(
        #[case] bond: ExtendedBond,
        #[case] a: u32,
        #[case] b: u32,
        #[case] expected: ExtendedBond,
    ) {
        assert_eq!(bond.update_atoms(a, b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::plain(ExtendedBond::new(0, 4, BondOrder::Single))]
    #[case::donating(ExtendedBond::new_dative(0, 4, BondOrder::Single, BondDonation::Donating))]
    #[case::narrowing_wedge(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }), ..ExtendedBond::new(0, 4, BondOrder::Single) })]
    fn test_extended_bond_update_atoms_identity(#[case] bond: ExtendedBond) {
        assert_eq!(bond.update_atoms(bond.start_atom(), bond.end_atom()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::none(ExtendedBond::new(0, 4, BondOrder::Single), None)]
    #[case::widening(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..ExtendedBond::new(0, 4, BondOrder::Single) }, Some(0))]
    #[case::narrowing(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Narrowing }), ..ExtendedBond::new(0, 4, BondOrder::Single) }, Some(4))]
    fn test_extended_bond_narrow_endpoint(#[case] bond: ExtendedBond, #[case] expected: Option<u32>) {
        assert_eq!(bond.narrow_endpoint(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::none(ExtendedBond::new(0, 4, BondOrder::Single), None)]
    #[case::widening(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Widening }), ..ExtendedBond::new(0, 4, BondOrder::Single) }, Some(4))]
    #[case::narrowing(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Up, taper: BondTaper::Narrowing }), ..ExtendedBond::new(0, 4, BondOrder::Single) }, Some(0))]
    fn test_extended_bond_wide_endpoint(#[case] bond: ExtendedBond, #[case] expected: Option<u32>) {
        assert_eq!(bond.wide_endpoint(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::normal(Bond::new(0, 1, BondOrder::Single), ExtendedBond::new(0, 1, BondOrder::Single))]
    #[case::wedge(Bond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }), ..Bond::new(0, 1, BondOrder::Single) }, ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }), ..ExtendedBond::new(0, 1, BondOrder::Single) })]
    fn test_bond_into_extended_bond(#[case] bond: Bond, #[case] expected: ExtendedBond) {
        let extended: ExtendedBond = bond.into();
        assert_eq!(extended, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::normal(ExtendedBond::new(0, 1, BondOrder::Double), Bond::new(0, 1, BondOrder::Double))]
    #[case::wedge(ExtendedBond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }), ..ExtendedBond::new(0, 1, BondOrder::Single) }, Bond { wedge: Some(BondWedge { orientation: BondOrientation::Down, taper: BondTaper::Narrowing }), ..Bond::new(0, 1, BondOrder::Single) })]
    fn test_extended_bond_try_into_bond(#[case] extended: ExtendedBond, #[case] expected: Bond) {
        let bond: Bond = extended.try_into().unwrap();
        assert_eq!(bond, expected);
    }

    #[rstest]
    #[case::query(ExtendedBond::new(0, 1, BondOrder::Any))]
    #[case::topology(ExtendedBond { atoms: AtomPair::new(0, 1), order: BondOrder::Single, topology: Some(BondTopology::Ring), donation: None, noncovalent: None,
                                    charge: None, unpaired_electrons: None, multiplicity: None, stereo: None, direction: None, wedge: None, ring: None,
                                    reacting_center: None, properties: HashMap::new(), span: None })]
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
    #[case::topology(ExtendedBond { atoms: AtomPair::new(0, 1), order: BondOrder::Single, topology: Some(BondTopology::Ring), donation: None, noncovalent: None,
                                    charge: None, unpaired_electrons: None, multiplicity: None, stereo: None, direction: None, wedge: None, ring: None,
                                    reacting_center: None, properties: HashMap::new(), span: None }, true)]
    fn test_extended_bond_has_extended_features(
        #[case] extended: ExtendedBond,
        #[case] expected: bool,
    ) {
        assert_eq!(extended.has_extended_features(), expected);
    }
}
