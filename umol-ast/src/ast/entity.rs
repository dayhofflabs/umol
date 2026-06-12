//! Molecule entity types.

use super::ids::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// A typed reference to any entity in a molecule — the variant is the kind, the
/// payload its id. General-purpose (coloring, symmetry analysis, constraints).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
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
