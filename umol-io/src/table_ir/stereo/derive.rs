//! Bond frames from supplied geometry and explicit CTfile/CX annotations.

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
pub(crate) fn derive_stereo_bonds<B>(
    atom_count: usize,
    bonds: &[B],
    bond_fields: impl Fn(&B) -> (AtomPair, BondOrder, Option<BondWedge>),
    positions: Option<&[Point3D]>,
    bond_stereo_assertions: &[(u32, BondStereo)],
) -> Result<Vec<StereoBond>, StereoDerivationError> {
    let mut codes = vec![None; bonds.len()];
    for &(bond, code) in bond_stereo_assertions {
        let Some((_, order, _)) = bonds.get(bond as usize).map(&bond_fields) else {
            return Err(StereoDerivationError::BondIndexOutOfBounds { bond });
        };
        if order != BondOrder::Double {
            return Err(StereoDerivationError::UnsupportedSite { bond });
        }
        if codes[bond as usize].is_some_and(|previous| previous != code) {
            return Err(StereoDerivationError::ConflictingConfiguration { bond });
        }
        codes[bond as usize] = Some(code);
    }
    let neighbors = AtomNeighbors::new(atom_count, bonds.iter().map(|bond| bond_fields(bond).0));
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
        let mut partners = neighbors
            .neighbors(atom)
            .iter()
            .filter(|neighbor| bond_fields(&bonds[neighbor.bond as usize]).1 == BondOrder::Double);
        if let (Some(partner), None) = (partners.next(), partners.next()) {
            let bond = partner.bond;
            if codes[bond as usize].is_some_and(|code| code != BondStereo::Either) {
                return Err(StereoDerivationError::ConflictingConfiguration { bond });
            }
            codes[bond as usize] = Some(BondStereo::Either);
        }
    }
    let mut frames = Vec::new();
    for (bond, (atoms, order, _)) in bonds.iter().map(&bond_fields).enumerate() {
        if order != BondOrder::Double {
            continue;
        }
        let code = codes[bond];
        if code == Some(BondStereo::Either) {
            frames.push(StereoBond {
                bond: bond as u32,
                configuration: BondConfiguration::Either,
            });
            continue;
        }
        if code.is_none() && positions.is_none() {
            continue;
        }
        for atom in [atoms.first(), atoms.second()] {
            if atom as usize >= atom_count {
                return Err(StereoDerivationError::AtomIndexOutOfBounds { atom });
            }
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
        let cumulated = |endpoint| {
            neighbors.neighbors(endpoint).iter().any(|neighbor| {
                neighbor.bond as usize != bond
                    && bond_fields(&bonds[neighbor.bond as usize]).1 == BondOrder::Double
            })
        };
        if first.is_empty()
            || second.is_empty()
            || first.len() > 2
            || second.len() > 2
            || atoms.first() == atoms.second()
            || first.contains(&atoms.first())
            || second.contains(&atoms.second())
            || first.iter().any(|atom| second.contains(atom))
            || cumulated(atoms.first())
            || cumulated(atoms.second())
        {
            if code.is_some() {
                return Err(StereoDerivationError::UnsupportedSite { bond: bond as u32 });
            }
            continue;
        }
        let geometry = positions
            .map(|positions| relation_from_positions(atoms, &first, &second, positions))
            .transpose()?
            .flatten();
        let explicit = match code {
            Some(BondStereo::Cis) => Some(BondRelation::SameSide),
            Some(BondStereo::Trans) => Some(BondRelation::OppositeSide),
            _ => None,
        };
        if let (Some(explicit), Some(geometry)) = (explicit, geometry) {
            if explicit != geometry {
                return Err(StereoDerivationError::ConflictingConfiguration { bond: bond as u32 });
            }
        }
        if let Some(relation) = explicit.or(geometry) {
            frames.push(StereoBond {
                bond: bond as u32,
                configuration: BondConfiguration::Framed {
                    references: [first[0], second[0]],
                    relation,
                },
            });
        }
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

#[cfg(all(test, feature = "proptest"))]
mod properties;
