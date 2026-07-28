//! Molecule entity types.

use std::fmt;

use strum::{EnumCount, EnumDiscriminants, FromRepr};

use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// A typed reference to any entity in a molecule — the variant is the kind, the
/// payload its id. General-purpose (coloring, symmetry analysis, constraints).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, EnumDiscriminants)]
#[strum_discriminants(name(EntityKind))]
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

impl fmt::Display for Entity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Entity::Atom(id) => write!(f, "atom {id}"),
            Entity::Bond(id) => write!(f, "bond {id}"),
            Entity::DativeBond(id) => write!(f, "dative bond {id}"),
            Entity::AromaticSystem(id) => write!(f, "aromatic system {id}"),
            Entity::MulticenterBond(id) => write!(f, "multicenter bond {id}"),
            Entity::NoncovalentBond(id) => write!(f, "noncovalent bond {id}"),
            Entity::StereoAtom(id) => write!(f, "stereo atom {id}"),
            Entity::StereoBond(id) => write!(f, "stereo bond {id}"),
        }
    }
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

impl TryFrom<u8> for EntityKind {
    type Error = ();

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(EntityKind::Atom),
            1 => Ok(EntityKind::Bond),
            2 => Ok(EntityKind::DativeBond),
            3 => Ok(EntityKind::AromaticSystem),
            4 => Ok(EntityKind::MulticenterBond),
            5 => Ok(EntityKind::NoncovalentBond),
            6 => Ok(EntityKind::StereoAtom),
            7 => Ok(EntityKind::StereoBond),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(1)), "atom 1")]
    #[case::bond(Entity::Bond(BondId(2)), "bond 2")]
    #[case::dative_bond(Entity::DativeBond(DativeBondId(3)), "dative bond 3")]
    #[case::aromatic_system(Entity::AromaticSystem(AromaticSystemId(4)), "aromatic system 4")]
    #[case::multicenter_bond(Entity::MulticenterBond(MulticenterBondId(5)), "multicenter bond 5")]
    #[case::noncovalent_bond(Entity::NoncovalentBond(NoncovalentBondId(6)), "noncovalent bond 6")]
    #[case::stereo_atom(Entity::StereoAtom(StereoAtomId(7)), "stereo atom 7")]
    #[case::stereo_bond(Entity::StereoBond(StereoBondId(8)), "stereo bond 8")]
    fn test_entity_display(#[case] entity: Entity, #[case] expected: &str) {
        assert_eq!(entity.to_string(), expected);
    }

    #[rstest]
    #[case::variant_and_id(
        vec![
            Entity::StereoBond(StereoBondId(0)),
            Entity::Atom(AtomId(2)),
            Entity::StereoAtom(StereoAtomId(0)),
            Entity::NoncovalentBond(NoncovalentBondId(0)),
            Entity::MulticenterBond(MulticenterBondId(0)),
            Entity::AromaticSystem(AromaticSystemId(0)),
            Entity::DativeBond(DativeBondId(0)),
            Entity::Bond(BondId(0)),
            Entity::Atom(AtomId(1)),
        ],
        vec![
            Entity::Atom(AtomId(1)),
            Entity::Atom(AtomId(2)),
            Entity::Bond(BondId(0)),
            Entity::DativeBond(DativeBondId(0)),
            Entity::AromaticSystem(AromaticSystemId(0)),
            Entity::MulticenterBond(MulticenterBondId(0)),
            Entity::NoncovalentBond(NoncovalentBondId(0)),
            Entity::StereoAtom(StereoAtomId(0)),
            Entity::StereoBond(StereoBondId(0)),
        ],
    )]
    fn test_entity_cmp(#[case] mut input: Vec<Entity>, #[case] expected: Vec<Entity>) {
        input.sort();

        assert_eq!(input, expected);
    }
}
