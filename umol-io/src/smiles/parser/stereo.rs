//! Double-bond frames derived from completed lexical direction markers.

use std::cell::{Cell, OnceCell};

use smallvec::SmallVec;

use super::super::error::ParseError;
use crate::table_ir::{
    AtomNeighbors, AtomPair, BondConfiguration, BondDirection, BondOrder, BondRelation, StereoBond,
};

pub(super) struct DirectionMarker {
    pub(super) direction: BondDirection,
    participating: Cell<bool>,
}

impl DirectionMarker {
    pub(super) fn new(direction: BondDirection) -> Self {
        Self {
            direction,
            participating: Cell::new(false),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub(super) enum DirectionError {
    AtomIndexOutOfBounds { atom: u32 },
    DanglingBondDirection { bond: u32 },
    CisTransConflict { atom: u32 },
    UnsupportedSite { bond: u32 },
}

/// Derive local frames in bond-table order. Directions use the first endpoint's viewpoint.
/// Accepts fresh markers and marks participation without consuming their directions.
/// Reuse the lookup only with the same atom count, bond-table indices, and endpoints.
pub(super) fn derive_stereo_bonds<B>(
    atom_count: usize,
    bonds: &[B],
    bond_fields: impl Fn(&B) -> (AtomPair, BondOrder, Option<&DirectionMarker>),
    neighbors: &OnceCell<AtomNeighbors>,
) -> Result<Vec<StereoBond>, DirectionError> {
    if !bonds.iter().any(|bond| bond_fields(bond).2.is_some()) {
        return Ok(Vec::new());
    }
    let mut frames = Vec::new();
    for (bond, (atoms, order, _)) in bonds.iter().map(&bond_fields).enumerate() {
        for atom in [atoms.first(), atoms.second()] {
            if atom as usize >= atom_count {
                return Err(DirectionError::AtomIndexOutOfBounds { atom });
            }
        }
        if order != BondOrder::Double {
            continue;
        }
        let neighbors = neighbors.get_or_init(|| {
            AtomNeighbors::new(atom_count, bonds.iter().map(|bond| bond_fields(bond).0))
        });
        let marked = |endpoint, other| {
            neighbors.neighbors(endpoint).iter().any(|neighbor| {
                neighbor.atom != other
                    && bond_fields(&bonds[neighbor.bond as usize]).1 == BondOrder::Single
                    && bond_fields(&bonds[neighbor.bond as usize]).2.is_some()
            })
        };
        if !marked(atoms.first(), atoms.second()) && !marked(atoms.second(), atoms.first()) {
            continue;
        }
        let substituents = |endpoint, other| {
            let mut atoms = SmallVec::<[u32; 2]>::new();
            let mut excess = false;
            for neighbor in neighbors.neighbors(endpoint) {
                if neighbor.atom == other || atoms.contains(&neighbor.atom) {
                    continue;
                }
                if atoms.len() == 2 {
                    excess = true;
                    break;
                }
                atoms.push(neighbor.atom);
            }
            atoms.sort_unstable();
            (atoms, excess)
        };
        let (first, first_excess) = substituents(atoms.first(), atoms.second());
        let (second, second_excess) = substituents(atoms.second(), atoms.first());
        if first.is_empty() || second.is_empty() {
            continue;
        }
        let cumulated = |endpoint| {
            neighbors.neighbors(endpoint).iter().any(|neighbor| {
                neighbor.bond as usize != bond
                    && bond_fields(&bonds[neighbor.bond as usize]).1 == BondOrder::Double
            })
        };
        if atoms.first() == atoms.second()
            || first.contains(&atoms.first())
            || second.contains(&atoms.second())
            || first_excess
            || second_excess
            || first.iter().any(|atom| second.contains(atom))
            || cumulated(atoms.first())
            || cumulated(atoms.second())
        {
            return Err(DirectionError::UnsupportedSite { bond: bond as u32 });
        }
        let side = |endpoint, other, reference| {
            let mut direction = None;
            for neighbor in neighbors.neighbors(endpoint) {
                if neighbor.atom == other {
                    continue;
                }
                let (pair, order, marker) = bond_fields(&bonds[neighbor.bond as usize]);
                if order != BondOrder::Single {
                    continue;
                }
                let Some(marker) = marker else { continue };
                let mut marker_direction = marker.direction;
                if pair.first() != endpoint {
                    marker_direction = marker_direction.flip();
                }
                if neighbor.atom != reference {
                    marker_direction = marker_direction.flip();
                }
                if direction.is_some_and(|previous| previous != marker_direction) {
                    return Err(DirectionError::CisTransConflict { atom: endpoint });
                }
                direction = Some(marker_direction);
                marker.participating.set(true);
            }
            Ok(direction)
        };
        let first_direction = side(atoms.first(), atoms.second(), first[0])?;
        let second_direction = side(atoms.second(), atoms.first(), second[0])?;
        if let (Some(first_direction), Some(second_direction)) = (first_direction, second_direction)
        {
            frames.push(StereoBond {
                bond: bond as u32,
                configuration: BondConfiguration::Framed {
                    references: [first[0], second[0]],
                    relation: if first_direction == second_direction {
                        BondRelation::SameSide
                    } else {
                        BondRelation::OppositeSide
                    },
                },
            });
        }
    }
    for (bond, (_, _, direction)) in bonds.iter().map(&bond_fields).enumerate() {
        if direction.is_some_and(|marker| !marker.participating.get()) {
            return Err(DirectionError::DanglingBondDirection { bond: bond as u32 });
        }
    }
    Ok(frames)
}

#[cfg(test)]
mod tests;

#[cfg(all(test, feature = "proptest"))]
mod properties;

impl From<DirectionError> for ParseError {
    fn from(error: DirectionError) -> Self {
        match error {
            DirectionError::AtomIndexOutOfBounds { atom } => {
                Self::AtomIndexOutOfBounds { atom_idx: atom }
            }
            DirectionError::DanglingBondDirection { bond } => Self::DanglingBondDirection { bond },
            DirectionError::CisTransConflict { atom } => Self::CisTransConflict { atom },
            DirectionError::UnsupportedSite { bond } => Self::UnsupportedStereoBond { bond },
        }
    }
}
