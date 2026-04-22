//! Per-atom constraints.

use strum::{EnumCount, EnumDiscriminants};

use crate::ast::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(AtomConstraintKind), derive(Hash, EnumCount))]
#[repr(u8)]
pub enum AtomConstraint {
    Valence(ValueAst),
    AromaticValence(AromaticValenceConstraint),
    MulticenterValence(MulticenterValenceConstraint),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    Degree(ValueAst),
    Connectivity(ValueAst),
    RingConnectivity(ValueAst),
    TotalHydrogens(ValueAst),
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl AtomConstraint {
    pub fn kind(&self) -> AtomConstraintKind {
        self.into()
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum AromaticValenceConstraint {
    #[default]
    Undetermined,
    NotAromatic,
    Aromatic(ValueAst),
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub enum MulticenterValenceConstraint {
    #[default]
    Undetermined,
    NotMulticenter,
    Multicenter(ValueAst),
}

/// Per-atom constraint slotmap. Fixed-size array indexed by
/// [`AtomConstraintKind`]; each slot holds at most one constraint of that
/// kind. O(1) access and update by kind.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AtomConstraints {
    slots: [Option<AtomConstraint>; AtomConstraintKind::COUNT],
}

impl Default for AtomConstraints {
    fn default() -> Self {
        Self {
            slots: [const { None }; AtomConstraintKind::COUNT],
        }
    }
}

impl AtomConstraints {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn is_empty(&self) -> bool {
        self.slots.iter().all(Option::is_none)
    }

    pub fn len(&self) -> usize {
        self.slots.iter().filter(|s| s.is_some()).count()
    }

    pub fn contains(&self, kind: AtomConstraintKind) -> bool {
        self.slots[kind as usize].is_some()
    }

    pub fn get(&self, kind: AtomConstraintKind) -> Option<&AtomConstraint> {
        self.slots[kind as usize].as_ref()
    }

    pub fn get_mut(&mut self, kind: AtomConstraintKind) -> Option<&mut AtomConstraint> {
        self.slots[kind as usize].as_mut()
    }

    /// Insert a constraint in its kind's slot, returning the previous occupant.
    pub fn set(&mut self, constraint: AtomConstraint) -> Option<AtomConstraint> {
        let slot = &mut self.slots[constraint.kind() as usize];
        slot.replace(constraint)
    }

    pub fn remove(&mut self, kind: AtomConstraintKind) -> Option<AtomConstraint> {
        self.slots[kind as usize].take()
    }

    pub fn iter(&self) -> impl Iterator<Item = &AtomConstraint> {
        self.slots.iter().filter_map(Option::as_ref)
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut AtomConstraint> {
        self.slots.iter_mut().filter_map(Option::as_mut)
    }
}

impl FromIterator<AtomConstraint> for AtomConstraints {
    fn from_iter<I: IntoIterator<Item = AtomConstraint>>(iter: I) -> Self {
        let mut out = Self::new();
        for c in iter {
            out.set(c);
        }
        out
    }
}
