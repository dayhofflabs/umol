//! Molecule entity types.

use strum::{EnumCount, EnumDiscriminants, FromRepr};

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// A typed reference to any entity in a molecule — the variant is the kind, the
/// payload its id. General-purpose (coloring, symmetry analysis, constraints).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, EnumDiscriminants)]
#[strum_discriminants(name(EntityKind), repr(u8))]
#[strum_discriminants(derive(Hash, EnumCount, FromRepr))]
pub enum Entity {
    Atom(AtomId),
    Bond(BondId),
    DativeBond(DativeBondId),
    AromaticSystem(AromaticSystemId),
    MulticenterBond(MulticenterBondId),
    NoncovalentBond(NoncovalentBondId),
    StereoAtom(StereoAtomId),
    StereoBond(StereoBondId),
}

impl Entity {
    pub fn kind(self) -> EntityKind {
        EntityKind::from(self)
    }

    /// Numeric index within entity kind.
    pub fn id_index(self) -> usize {
        match self {
            Entity::Atom(i) => i.index(),
            Entity::Bond(i) => i.index(),
            Entity::DativeBond(i) => i.index(),
            Entity::AromaticSystem(i) => i.index(),
            Entity::MulticenterBond(i) => i.index(),
            Entity::NoncovalentBond(i) => i.index(),
            Entity::StereoAtom(i) => i.index(),
            Entity::StereoBond(i) => i.index(),
        }
    }
}

impl EntityKind {
    /// Reconstruct entity of this kind with the given id index.
    pub fn with_id(self, id: u32) -> Entity {
        match self {
            EntityKind::Atom => Entity::Atom(AtomId(id)),
            EntityKind::Bond => Entity::Bond(BondId(id)),
            EntityKind::DativeBond => Entity::DativeBond(DativeBondId(id)),
            EntityKind::AromaticSystem => Entity::AromaticSystem(AromaticSystemId(id)),
            EntityKind::MulticenterBond => Entity::MulticenterBond(MulticenterBondId(id)),
            EntityKind::NoncovalentBond => Entity::NoncovalentBond(NoncovalentBondId(id)),
            EntityKind::StereoAtom => Entity::StereoAtom(StereoAtomId(id)),
            EntityKind::StereoBond => Entity::StereoBond(StereoBondId(id)),
        }
    }
}
