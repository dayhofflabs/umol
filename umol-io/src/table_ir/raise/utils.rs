//! Raise-time stereo helpers for `super`: ligand orderings, directional-bond faces, and the
//! capability/validation predicates the `raise_*` functions build on.

use umol_geometric_core::{
    complementary_direction, same_side_of_axis, signed_volume, Point3D, AXIS_SIDE_TOLERANCE,
};
use umol_graph_ir::ir::NoncovalentBondKind;

use super::RaiseError;
use crate::table_ir::bond::{BondNoncovalent as TableNoncovalent, BondOrder as TableBondOrder};
use crate::table_ir::{AtomNeighbors, BondDirection, BondOrientation, Molecule as TableMolecule};

pub(super) fn noncovalent_kind(kind: TableNoncovalent) -> NoncovalentBondKind {
    match kind {
        TableNoncovalent::Hydrogen => NoncovalentBondKind::HydrogenBond,
    }
}

/// Ligand stereo ordering for raise operation (atom or virtual ligand). Virtual
/// ligand does not distinguish between implicit H or lone pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum StereoLigand {
    Atom(usize),
    Virtual(usize),
}

/// Halfplane of the plane of the double bond, split by the bond axis.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum StereoHalfplane {
    Top,
    Bottom,
}

impl StereoHalfplane {
    fn flip(self) -> Self {
        match self {
            Self::Top => Self::Bottom,
            Self::Bottom => Self::Top,
        }
    }
}

/// Out-of-plane direction.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum StereoOutofPlane {
    Front,
    Back,
}

/// Atom adjacent to stereogenic double bond. Second virtual ligand is
/// added if atom has one substituent. Halfplane of the second ligand is
/// flipped from the first.
pub(super) struct StereoBondAtom {
    pub(super) first_ligand: StereoLigand,
    pub(super) second_ligand: StereoLigand,
    pub(super) first_halfplane: StereoHalfplane,
}

/// Neighbor atom ordering of `atom_idx`, by ascending atom index.
fn atom_ordering(neighbors: &AtomNeighbors, atom_idx: usize) -> Vec<usize> {
    let mut indices: Vec<usize> = neighbors
        .neighbors(atom_idx as u32)
        .iter()
        .map(|neighbor| neighbor.atom as usize)
        .collect();
    indices.sort_unstable();
    indices.dedup();
    indices
}

/// Number of distinct atoms neighboring `atom_idx`.
pub(super) fn neighbor_count(neighbors: &AtomNeighbors, atom_idx: usize) -> usize {
    neighbors.degree(atom_idx as u32)
}

/// Neighbor atom ordering of `atom_idx` by bond ordering (used by SMILES, which refers to it as parse ordering).
/// Neighbor atoms appear in the order of incident bonds (including ring-closure indices).
/// Can differ from the ascending atom index order when rings are present.
fn bond_neighbor_ordering(neighbors: &AtomNeighbors, atom_idx: usize) -> Vec<usize> {
    let mut indices = Vec::new();
    for other in neighbors
        .neighbors(atom_idx as u32)
        .iter()
        .map(|neighbor| neighbor.atom as usize)
    {
        if !indices.contains(&other) {
            indices.push(other);
        }
    }
    indices
}

/// Ligand ordering used in tetrahedral stereo constraints (#T): neighbors ascending as `Atom`,
/// then at most one `Virtual`. More than one virtual ligand is disallowed by `validate_tetrahedral_geometry`.
pub(super) fn tetrahedral_ligand_ordering(
    neighbors: &AtomNeighbors,
    atom_idx: usize,
) -> Vec<StereoLigand> {
    let mut ordering: Vec<StereoLigand> = atom_ordering(neighbors, atom_idx)
        .into_iter()
        .map(StereoLigand::Atom)
        .collect();
    if ordering.len() == 3 {
        ordering.push(StereoLigand::Virtual(atom_idx));
    }
    ordering
}

/// SMILES/SMARTS tetrahedral ligand ordering, FirstNeighborToward: neighbors in parse order, virtual
/// ligand is first if `atom_idx` opens the SMILES, else second.
pub(super) fn first_neighbor_toward_ordering(
    neighbors: &AtomNeighbors,
    atom_idx: usize,
) -> Vec<StereoLigand> {
    let mut ordering: Vec<StereoLigand> = bond_neighbor_ordering(neighbors, atom_idx)
        .into_iter()
        .map(StereoLigand::Atom)
        .collect();
    if ordering.len() == 3 {
        let ligand_idx = if atom_idx > 0 { 1 } else { 0 };
        ordering.insert(ligand_idx, StereoLigand::Virtual(atom_idx));
    }
    ordering
}

/// Tetrahedral stereo coset index from wedge bonds at `atom_idx`, or `None` when the projected
/// ligand positions are all at the center or collinear, so that no winding exists.
pub(super) fn coset_from_wedge_winding(
    atom_idx: usize,
    ordering: &[StereoLigand],
    wedged: usize,
    positions: &[Point3D],
    outofplane: StereoOutofPlane,
) -> Option<usize> {
    let z = if outofplane == StereoOutofPlane::Front {
        1.0
    } else {
        -1.0
    };
    let center_position = positions[atom_idx];
    let neighbor_positions: Vec<Point3D> = ordering
        .iter()
        .filter_map(|&ligand| match ligand {
            StereoLigand::Atom(index) => Some(positions[index]),
            StereoLigand::Virtual(_) => None,
        })
        .collect();
    let virtual_position = complementary_direction(center_position, &neighbor_positions);
    let points: Vec<Point3D> = ordering
        .iter()
        .map(|&ligand| match ligand {
            StereoLigand::Atom(index) if index == wedged => {
                Point3D::new(positions[wedged].x, positions[wedged].y, z)
            }
            StereoLigand::Atom(index) => Point3D::new(positions[index].x, positions[index].y, 0.0),
            StereoLigand::Virtual(_) => virtual_position,
        })
        .collect();
    let scale = neighbor_positions
        .iter()
        .map(|position| {
            let (dx, dy) = (
                position.x - center_position.x,
                position.y - center_position.y,
            );
            (dx * dx + dy * dy).sqrt()
        })
        .fold(0.0, f64::max);
    let volume = signed_volume(points[0], points[1], points[2], points[3]);
    if scale == 0.0 || volume.abs() <= AXIS_SIDE_TOLERANCE * scale * scale {
        return None;
    }
    // umol convention (matching the SMILES `@` = anticlockwise = coset 0 path): the ascending-index
    // ligands of coset 0 have a negative signed volume.
    Some(if volume < 0.0 { 0 } else { 1 })
}

/// Wide endpoints of the definite wedges whose narrow end is `atom_idx`, each with its
/// out-of-plane direction. A wedge describes only its narrow endpoint; `Either` wedges are read
/// by `has_either_wedge`.
pub(super) fn wedge_bond_neighbors(
    mol: &TableMolecule,
    neighbors: &AtomNeighbors,
    atom_idx: usize,
) -> Vec<(usize, StereoOutofPlane)> {
    neighbors
        .neighbors(atom_idx as u32)
        .iter()
        .filter_map(|neighbor| {
            let bond = &mol.bonds[neighbor.bond as usize];
            if bond.narrow_endpoint() != Some(atom_idx as u32) {
                return None;
            }
            match bond.wedge?.orientation {
                BondOrientation::Up => Some((neighbor.atom as usize, StereoOutofPlane::Front)),
                BondOrientation::Down => Some((neighbor.atom as usize, StereoOutofPlane::Back)),
                BondOrientation::Either
                | BondOrientation::EitherUp
                | BondOrientation::EitherDown => None,
            }
        })
        .collect()
}

/// The partner of `atom_idx`'s double bond when it has exactly one.
pub(super) fn double_bond_partner(
    mol: &TableMolecule,
    neighbors: &AtomNeighbors,
    atom_idx: usize,
) -> Option<usize> {
    let mut partners = neighbors
        .neighbors(atom_idx as u32)
        .iter()
        .filter(|neighbor| mol.bonds[neighbor.bond as usize].order == TableBondOrder::Double);
    match (partners.next(), partners.next()) {
        (Some(partner), None) => Some(partner.atom as usize),
        _ => None,
    }
}

/// Whether an `Either` wedge (MOL code 4, CXSMILES `w:`, `wU:`, `wD:`) has its narrow end at
/// `atom_idx`. At an atom with exactly one double bond the mark is the drawing convention for an
/// unknown configuration of that double bond; elsewhere it asserts a stereo center of unknown
/// configuration.
pub(super) fn has_either_wedge(
    mol: &TableMolecule,
    neighbors: &AtomNeighbors,
    atom_idx: usize,
) -> bool {
    neighbors.neighbors(atom_idx as u32).iter().any(|neighbor| {
        let bond = &mol.bonds[neighbor.bond as usize];
        bond.narrow_endpoint() == Some(atom_idx as u32)
            && matches!(
                bond.wedge.map(|wedge| wedge.orientation),
                Some(
                    BondOrientation::Either
                        | BondOrientation::EitherUp
                        | BondOrientation::EitherDown
                )
            )
    })
}

/// Validate that tetrahedral stereo has 3 or 4 neighbors.
pub(super) fn validate_tetrahedral_geometry(
    neighbors: &AtomNeighbors,
    atom_idx: usize,
) -> Result<(), RaiseError> {
    let count = neighbor_count(neighbors, atom_idx);
    if count == 3 || count == 4 {
        Ok(())
    } else {
        Err(RaiseError::TetrahedralLigandCount {
            atom: atom_idx,
            count,
        })
    }
}

/// Validate that a directional bond (`/`,`\`) is adjacent to a cis/trans-capable double bond.
/// Returns `Ok(())` for any non-directional bond.
pub(super) fn validate_bond_direction(
    mol: &TableMolecule,
    neighbors: &AtomNeighbors,
    bond_idx: usize,
) -> Result<(), RaiseError> {
    let bond = &mol.bonds[bond_idx];
    if bond.order != TableBondOrder::Single || bond.direction.is_none() {
        return Ok(());
    }
    let flanks_capable = [bond.start_atom(), bond.end_atom()]
        .into_iter()
        .any(|atom| {
            neighbors.neighbors(atom).iter().any(|neighbor| {
                mol.bonds[neighbor.bond as usize].order == TableBondOrder::Double
                    && cis_trans_capable(neighbors, atom as usize, neighbor.atom as usize)
            })
        });
    if flanks_capable {
        Ok(())
    } else {
        Err(RaiseError::DanglingBondDirection { bond: bond_idx })
    }
}

/// Double bond is cis-trans capable iff both ends have distinct substituents.
pub(super) fn cis_trans_capable(neighbors: &AtomNeighbors, atom_1: usize, atom_2: usize) -> bool {
    let side_1: Vec<_> = atom_ordering(neighbors, atom_1)
        .into_iter()
        .filter(|&atom| atom != atom_2)
        .collect();
    let side_2: Vec<_> = atom_ordering(neighbors, atom_2)
        .into_iter()
        .filter(|&atom| atom != atom_1)
        .collect();
    !side_1.is_empty() && !side_2.is_empty() && side_1.iter().all(|atom| !side_2.contains(atom))
}

/// Arrangement of the bond atom `atom_idx` of stereogenic double bond. Errors when its markers disagree.
pub(super) fn cis_trans_side(
    mol: &TableMolecule,
    neighbors: &AtomNeighbors,
    atom_idx: usize,
    other_atom_idx: usize,
) -> Result<Option<StereoBondAtom>, RaiseError> {
    let substituents: Vec<usize> = atom_ordering(neighbors, atom_idx)
        .into_iter()
        .filter(|&n| n != other_atom_idx)
        .collect();
    let Some(&first) = substituents.first() else {
        return Ok(None);
    };
    let first_ligand = StereoLigand::Atom(first);
    let second_ligand = substituents
        .get(1)
        .map_or(StereoLigand::Virtual(atom_idx), |&second| {
            StereoLigand::Atom(second)
        });
    // The first ligand's face: from the bond toward it, or the flipped bond toward the geminal second.
    let toward_first = direction(mol, neighbors, atom_idx, first);
    let toward_second = substituents
        .get(1)
        .and_then(|&second| direction(mol, neighbors, atom_idx, second))
        .map(StereoHalfplane::flip);
    let first_halfplane = match (toward_first, toward_second) {
        (Some(a), Some(b)) if a != b => {
            return Err(RaiseError::CisTransConflict { atom: atom_idx })
        }
        (Some(halfplane), _) | (_, Some(halfplane)) => halfplane,
        (None, None) => return Ok(None),
    };
    Ok(Some(StereoBondAtom {
        first_ligand,
        second_ligand,
        first_halfplane,
    }))
}

/// Arrangements of both atoms of the double bond read from coordinates, as the specification
/// prescribes for stereo code 0. `None` when either atom has no substituent or is a cumulated
/// center, whose only other bond is a second double bond and which therefore has no plane of
/// substituents, and whenever the drawing does not settle the configuration: coincident bond
/// atoms, a substituent on the bond axis, or both substituents of one atom on one side of it. The
/// coordinates are supplementary and assert nothing they do not show.
pub(super) fn cis_trans_sides_from_positions(
    mol: &TableMolecule,
    neighbors: &AtomNeighbors,
    atom_1_idx: usize,
    atom_2_idx: usize,
    positions: &[Point3D],
) -> Option<(StereoBondAtom, StereoBondAtom)> {
    let substituents = |atom_idx: usize, other_atom_idx: usize| {
        let mut substituents: Vec<(usize, TableBondOrder)> = neighbors
            .neighbors(atom_idx as u32)
            .iter()
            .filter(|neighbor| neighbor.atom as usize != other_atom_idx)
            .map(|neighbor| {
                (
                    neighbor.atom as usize,
                    mol.bonds[neighbor.bond as usize].order,
                )
            })
            .collect();
        substituents.sort_unstable_by_key(|&(neighbor, _)| neighbor);
        substituents.dedup_by_key(|&mut (neighbor, _)| neighbor);
        substituents
    };
    let substituents_1 = substituents(atom_1_idx, atom_2_idx);
    let substituents_2 = substituents(atom_2_idx, atom_1_idx);
    let cumulated = |substituents: &[(usize, TableBondOrder)]| {
        matches!(substituents, [(_, TableBondOrder::Double)])
    };
    if substituents_1.is_empty()
        || substituents_2.is_empty()
        || cumulated(&substituents_1)
        || cumulated(&substituents_2)
    {
        return None;
    }
    let same_side = |first: usize, second: usize| {
        same_side_of_axis(
            positions[atom_1_idx],
            positions[atom_2_idx],
            positions[first],
            positions[second],
        )
    };
    // Two substituents of one atom must lie on opposite sides of the axis.
    for substituents in [&substituents_1, &substituents_2] {
        if let &[(first, _), (second, _)] = substituents.as_slice() {
            if same_side(first, second)? {
                return None;
            }
        }
    }
    let cis = same_side(substituents_1[0].0, substituents_2[0].0)?;
    let arrangement = |atom_idx: usize,
                       substituents: &[(usize, TableBondOrder)],
                       first_halfplane: StereoHalfplane| StereoBondAtom {
        first_ligand: StereoLigand::Atom(substituents[0].0),
        second_ligand: match substituents.get(1) {
            Some(&(second, _)) => StereoLigand::Atom(second),
            None => StereoLigand::Virtual(atom_idx),
        },
        first_halfplane,
    };
    Some((
        arrangement(atom_1_idx, &substituents_1, StereoHalfplane::Top),
        arrangement(
            atom_2_idx,
            &substituents_2,
            if cis {
                StereoHalfplane::Top
            } else {
                StereoHalfplane::Bottom
            },
        ),
    ))
}

/// Halfplane (top/bottom) of `other_atom_idx` viewed from `atom_idx`.
fn direction(
    mol: &TableMolecule,
    neighbors: &AtomNeighbors,
    atom_idx: usize,
    other_atom_idx: usize,
) -> Option<StereoHalfplane> {
    neighbors
        .neighbors(atom_idx as u32)
        .iter()
        .find_map(|neighbor| {
            if neighbor.atom as usize != other_atom_idx {
                return None;
            }
            let bond = &mol.bonds[neighbor.bond as usize];
            if bond.order != TableBondOrder::Single {
                return None;
            }
            let face = match bond.direction? {
                BondDirection::Rising => StereoHalfplane::Top,
                BondDirection::Falling => StereoHalfplane::Bottom,
            };
            Some(if bond.start_atom() as usize == atom_idx {
                face
            } else {
                face.flip()
            })
        })
}
