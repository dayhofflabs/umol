//! Double-bond frames derived from completed lexical direction markers.

use super::super::error::ParseError;
use crate::table_ir::{
    AtomNeighbors, AtomPair, BondConfiguration, BondDirection, BondOrder, BondRelation, StereoBond,
};

#[derive(Debug, PartialEq, Eq)]
pub(super) enum DirectionError {
    AtomIndexOutOfBounds { atom: u32 },
    DanglingBondDirection { bond: u32 },
    CisTransConflict { atom: u32 },
    UnsupportedSite { bond: u32 },
}

/// Derive local frames in bond-table order. Directions use the first endpoint's viewpoint.
pub(super) fn derive_stereo_bonds(
    atom_count: usize,
    bonds: &[(AtomPair, BondOrder, Option<BondDirection>)],
) -> Result<Vec<StereoBond>, DirectionError> {
    if !bonds.iter().any(|(_, _, direction)| direction.is_some()) {
        return Ok(Vec::new());
    }
    let neighbors = AtomNeighbors::new(atom_count, bonds.iter().map(|(atoms, _, _)| *atoms));
    let mut participating = vec![false; bonds.len()];
    let mut frames = Vec::new();
    for (bond, &(atoms, order, _)) in bonds.iter().enumerate() {
        for atom in [atoms.first(), atoms.second()] {
            if atom as usize >= atom_count {
                return Err(DirectionError::AtomIndexOutOfBounds { atom });
            }
        }
        if order != BondOrder::Double {
            continue;
        }
        let substituents = |endpoint, other| {
            let mut atoms: Vec<_> = neighbors
                .neighbors(endpoint)
                .iter()
                .filter(|neighbor| neighbor.atom != other)
                .map(|neighbor| neighbor.atom)
                .collect();
            atoms.sort_unstable();
            atoms.dedup();
            atoms
        };
        let first = substituents(atoms.first(), atoms.second());
        let second = substituents(atoms.second(), atoms.first());
        let marked = |endpoint, other| {
            neighbors.neighbors(endpoint).iter().any(|neighbor| {
                neighbor.atom != other
                    && bonds[neighbor.bond as usize].1 == BondOrder::Single
                    && bonds[neighbor.bond as usize].2.is_some()
            })
        };
        if !marked(atoms.first(), atoms.second()) && !marked(atoms.second(), atoms.first()) {
            continue;
        }
        if first.is_empty() || second.is_empty() {
            continue;
        }
        let cumulated = |endpoint| {
            neighbors.neighbors(endpoint).iter().any(|neighbor| {
                neighbor.bond as usize != bond
                    && bonds[neighbor.bond as usize].1 == BondOrder::Double
            })
        };
        if atoms.first() == atoms.second()
            || first.contains(&atoms.first())
            || second.contains(&atoms.second())
            || first.len() > 2
            || second.len() > 2
            || first.iter().any(|atom| second.contains(atom))
            || cumulated(atoms.first())
            || cumulated(atoms.second())
        {
            return Err(DirectionError::UnsupportedSite { bond: bond as u32 });
        }
        let mut side = |endpoint, other, reference| {
            let mut direction = None;
            for neighbor in neighbors.neighbors(endpoint) {
                if neighbor.atom == other {
                    continue;
                }
                let (pair, order, marker) = bonds[neighbor.bond as usize];
                if order != BondOrder::Single {
                    continue;
                }
                let Some(mut marker) = marker else { continue };
                if pair.first() != endpoint {
                    marker = marker.flip();
                }
                if neighbor.atom != reference {
                    marker = marker.flip();
                }
                if direction.is_some_and(|previous| previous != marker) {
                    return Err(DirectionError::CisTransConflict { atom: endpoint });
                }
                direction = Some(marker);
                participating[neighbor.bond as usize] = true;
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
    for (bond, (_, _, direction)) in bonds.iter().enumerate() {
        if direction.is_some() && !participating[bond] {
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
