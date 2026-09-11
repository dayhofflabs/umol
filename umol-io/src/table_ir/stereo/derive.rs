//! Bond frames from supplied geometry and explicit CTfile/CX annotations.

use std::cell::OnceCell;

use smallvec::SmallVec;
use umol_geometric_core::{same_side_of_axis, Point3D};

use crate::table_ir::{
    AtomNeighbors, AtomPair, BondConfiguration, BondOrder, BondOrientation, BondRelation,
    BondStereo, BondTaper, BondWedge, StereoBond,
};

#[derive(Debug, PartialEq, Eq)]
pub(crate) enum StereoDerivationError {
    BondIndexOutOfBounds { bond: u32 },
    AtomIndexOutOfBounds { atom: u32 },
    MissingPosition { atom: u32 },
    UnsupportedSite { bond: u32 },
    ConflictingConfiguration { bond: u32 },
}

/// Derive frames in table order, using supplied geometry only where no Either code applies.
/// Existing parser-produced frames are unique and ordered by bond; retain their references.
/// Reuse the lookup only with the same atom count, bond-table indices, and endpoints.
pub(crate) fn derive_stereo_bonds<B>(
    atom_count: usize,
    bonds: &[B],
    bond_fields: impl Fn(&B) -> (AtomPair, BondOrder, Option<BondWedge>),
    positions: Option<&[Point3D]>,
    mut frames: Vec<StereoBond>,
    mut bond_stereo_assertions: Vec<(u32, BondStereo)>,
    neighbors: &OnceCell<AtomNeighbors>,
) -> Result<Vec<StereoBond>, StereoDerivationError> {
    for &(bond, _) in &bond_stereo_assertions {
        let Some((_, order, _)) = bonds.get(bond as usize).map(&bond_fields) else {
            return Err(StereoDerivationError::BondIndexOutOfBounds { bond });
        };
        if order != BondOrder::Double {
            return Err(StereoDerivationError::UnsupportedSite { bond });
        }
    }
    for (atoms, _, wedge) in bonds.iter().map(&bond_fields) {
        let Some(wedge) = wedge.filter(|wedge| {
            matches!(
                wedge.orientation,
                BondOrientation::Either | BondOrientation::EitherUp | BondOrientation::EitherDown
            )
        }) else {
            continue;
        };
        let atom = match wedge.taper {
            BondTaper::Widening => atoms.first(),
            BondTaper::Narrowing => atoms.second(),
        };
        if atom as usize >= atom_count {
            return Err(StereoDerivationError::AtomIndexOutOfBounds { atom });
        }
        let neighbors = neighbors.get_or_init(|| {
            AtomNeighbors::new(atom_count, bonds.iter().map(|bond| bond_fields(bond).0))
        });
        let mut partners = neighbors
            .neighbors(atom)
            .iter()
            .filter(|neighbor| bond_fields(&bonds[neighbor.bond as usize]).1 == BondOrder::Double);
        if let (Some(partner), None) = (partners.next(), partners.next()) {
            let bond = partner.bond;
            bond_stereo_assertions.push((bond, BondStereo::Either));
        }
    }
    bond_stereo_assertions.sort_unstable_by_key(|&(bond, _)| bond);
    for pair in bond_stereo_assertions.windows(2) {
        if pair[0].0 == pair[1].0 && pair[0].1 != pair[1].1 {
            return Err(StereoDerivationError::ConflictingConfiguration { bond: pair[0].0 });
        }
    }
    bond_stereo_assertions.dedup();
    for frame in &frames {
        let Some((_, order, _)) = bonds.get(frame.bond as usize).map(&bond_fields) else {
            return Err(StereoDerivationError::BondIndexOutOfBounds { bond: frame.bond });
        };
        if order != BondOrder::Double {
            return Err(StereoDerivationError::UnsupportedSite { bond: frame.bond });
        }
    }
    let existing_count = frames.len();
    let mut frame_index = 0;
    let mut codes = bond_stereo_assertions.into_iter().peekable();
    for (bond, (atoms, order, _)) in bonds.iter().map(&bond_fields).enumerate() {
        if order != BondOrder::Double {
            continue;
        }
        let code = codes
            .next_if(|&(site, _)| site as usize == bond)
            .map(|(_, code)| code);
        let existing = if frame_index < existing_count && frames[frame_index].bond as usize == bond
        {
            let configuration = frames[frame_index].configuration;
            frame_index += 1;
            Some(configuration)
        } else {
            None
        };
        if code == Some(BondStereo::Either) || existing == Some(BondConfiguration::Either) {
            if code.is_some_and(|code| code != BondStereo::Either)
                || matches!(existing, Some(BondConfiguration::Framed { .. }))
            {
                return Err(StereoDerivationError::ConflictingConfiguration { bond: bond as u32 });
            }
            if existing.is_none() {
                frames.push(StereoBond {
                    bond: bond as u32,
                    configuration: BondConfiguration::Either,
                });
            }
            continue;
        }
        if code.is_none() && existing.is_none() && positions.is_none() {
            continue;
        }
        for atom in [atoms.first(), atoms.second()] {
            if atom as usize >= atom_count {
                return Err(StereoDerivationError::AtomIndexOutOfBounds { atom });
            }
        }
        let neighbors = neighbors.get_or_init(|| {
            AtomNeighbors::new(atom_count, bonds.iter().map(|bond| bond_fields(bond).0))
        });
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
        let cumulated = |endpoint| {
            neighbors.neighbors(endpoint).iter().any(|neighbor| {
                neighbor.bond as usize != bond
                    && bond_fields(&bonds[neighbor.bond as usize]).1 == BondOrder::Double
            })
        };
        if first.is_empty()
            || second.is_empty()
            || first_excess
            || second_excess
            || atoms.first() == atoms.second()
            || first.contains(&atoms.first())
            || second.contains(&atoms.second())
            || first.iter().any(|atom| second.contains(atom))
            || cumulated(atoms.first())
            || cumulated(atoms.second())
        {
            if code.is_some() || existing.is_some() {
                return Err(StereoDerivationError::UnsupportedSite { bond: bond as u32 });
            }
            continue;
        }
        let previous = if let Some(BondConfiguration::Framed {
            references,
            relation,
        }) = existing
        {
            if !first.contains(&references[0]) || !second.contains(&references[1]) {
                return Err(StereoDerivationError::UnsupportedSite { bond: bond as u32 });
            }
            Some(
                if (references[0] != first[0]) ^ (references[1] != second[0]) {
                    match relation {
                        BondRelation::SameSide => BondRelation::OppositeSide,
                        BondRelation::OppositeSide => BondRelation::SameSide,
                    }
                } else {
                    relation
                },
            )
        } else {
            None
        };
        let explicit = match code {
            Some(BondStereo::Cis) => Some(BondRelation::SameSide),
            Some(BondStereo::Trans) => Some(BondRelation::OppositeSide),
            _ => None,
        };
        if let (Some(explicit), Some(previous)) = (explicit, previous) {
            if explicit != previous {
                return Err(StereoDerivationError::ConflictingConfiguration { bond: bond as u32 });
            }
        }
        let explicit = explicit.or(previous);
        let geometry = positions
            .map(|positions| relation_from_positions(atoms, &first, &second, positions))
            .transpose()?
            .flatten();
        if let (Some(explicit), Some(geometry)) = (explicit, geometry) {
            if explicit != geometry {
                return Err(StereoDerivationError::ConflictingConfiguration { bond: bond as u32 });
            }
        }
        if let Some(relation) = explicit.or(geometry).filter(|_| existing.is_none()) {
            frames.push(StereoBond {
                bond: bond as u32,
                configuration: BondConfiguration::Framed {
                    references: [first[0], second[0]],
                    relation,
                },
            });
        }
    }
    if existing_count > 0 && frames.len() > existing_count {
        frames.sort_unstable_by_key(|frame| frame.bond);
    }
    Ok(frames)
}

fn relation_from_positions(
    atoms: AtomPair,
    first: &[u32],
    second: &[u32],
    positions: &[Point3D],
) -> Result<Option<BondRelation>, StereoDerivationError> {
    let same_side = |a, b| {
        let mut points = [Point3D::zero(); 4];
        for (point, atom) in points.iter_mut().zip([atoms.first(), atoms.second(), a, b]) {
            *point = *positions
                .get(atom as usize)
                .ok_or(StereoDerivationError::MissingPosition { atom })?;
        }
        if points
            .iter()
            .any(|p| !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite())
        {
            return Ok(None);
        }
        // Uniform scaling bounds intermediate products without changing the relative tolerance.
        let scale = points
            .iter()
            .flat_map(|p| [p.x.abs(), p.y.abs(), p.z.abs()])
            .fold(0.0, f64::max);
        if scale == 0.0 {
            return Ok(None);
        }
        let points = points.map(|p| Point3D::new(p.x / scale, p.y / scale, p.z / scale));
        Ok(same_side_of_axis(
            points[0], points[1], points[2], points[3],
        ))
    };
    for side in [first, second] {
        if let [a, b] = side {
            if same_side(*a, *b)? != Some(false) {
                return Ok(None);
            }
        }
    }
    Ok(same_side(first[0], second[0])?.map(|same| {
        if same {
            BondRelation::SameSide
        } else {
            BondRelation::OppositeSide
        }
    }))
}

#[cfg(test)]
mod tests;
