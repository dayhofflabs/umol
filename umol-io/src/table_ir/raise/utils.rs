//! Raise-time stereo helpers for `super`: ligand orderings, directional-bond faces, and the
//! capability/validation predicates the `raise_*` functions build on.

use umol_ast::ast::NoncovalentBondKind;
use umol_geometric_core::{complementary_direction, signed_volume, Point3D};

use super::RaiseError;
use crate::table_ir::bond::{BondNoncovalent as TableNoncovalent, BondOrder as TableBondOrder};
use crate::table_ir::{BondDirection, BondWedge, Molecule as TableMolecule};

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
fn atom_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<usize> {
    let mut indices: Vec<usize> = mol
        .bonds
        .iter()
        .filter_map(|bond| bond.atoms.other(atom_idx as u32))
        .map(|other| other as usize)
        .collect();
    indices.sort_unstable();
    indices.dedup();
    indices
}

/// Number of distinct atoms neighboring `atom_idx`.
pub(super) fn neighbor_count(mol: &TableMolecule, atom_idx: usize) -> usize {
    atom_ordering(mol, atom_idx).len()
}

/// Neighbor atom ordering of `atom_idx` by bond ordering (used by SMILES, which refers to it as parse ordering).
/// Neighbor atoms appear in the order of incident bonds (including ring-closure indices).
/// Can differ from the ascending atom index order when rings are present.
fn bond_neighbor_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<usize> {
    let mut indices = Vec::new();
    for other in mol
        .bonds
        .iter()
        .filter_map(|bond| bond.atoms.other(atom_idx as u32))
        .map(|other| other as usize)
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
    mol: &TableMolecule,
    atom_idx: usize,
) -> Vec<StereoLigand> {
    let mut ordering: Vec<StereoLigand> = atom_ordering(mol, atom_idx)
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
    mol: &TableMolecule,
    atom_idx: usize,
) -> Vec<StereoLigand> {
    let mut ordering: Vec<StereoLigand> = bond_neighbor_ordering(mol, atom_idx)
        .into_iter()
        .map(StereoLigand::Atom)
        .collect();
    if ordering.len() == 3 {
        let ligand_idx = if atom_idx > 0 { 1 } else { 0 };
        ordering.insert(ligand_idx, StereoLigand::Virtual(atom_idx));
    }
    ordering
}

/// MDL/CTFile tetradehral ligand ordering, LastNeighborAway: atom index ordering,
/// virtual ligand is last.
pub(super) fn last_neighbor_away_ordering(
    mol: &TableMolecule,
    atom_idx: usize,
) -> Vec<StereoLigand> {
    tetrahedral_ligand_ordering(mol, atom_idx)
}

/// Tetrahedral stereo coset index from wedge bonds at `atom_idx`.
pub(super) fn coset_from_wedge_winding(
    atom_idx: usize,
    ordering: &[StereoLigand],
    wedged: usize,
    positions: &[Point3D],
    outofplane: StereoOutofPlane,
) -> usize {
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
    // umol convention (matching the SMILES `@` = anticlockwise = coset 0 path): the ascending-index
    // ligands of coset 0 have a negative signed volume.
    if signed_volume(points[0], points[1], points[2], points[3]) < 0.0 {
        0
    } else {
        1
    }
}

/// Neighbors reached by wedge bonds at `atom_idx`, each with its out-of-plane direction.
pub(super) fn wedge_bond_neighbors(
    mol: &TableMolecule,
    atom_idx: usize,
) -> Vec<(usize, StereoOutofPlane)> {
    mol.bonds
        .iter()
        .filter_map(
            |bond| match (bond.atoms.other(atom_idx as u32), bond.wedge) {
                (Some(other), Some(BondWedge::Up)) => {
                    Some((other as usize, StereoOutofPlane::Front))
                }
                (Some(other), Some(BondWedge::Down)) => {
                    Some((other as usize, StereoOutofPlane::Back))
                }
                _ => None,
            },
        )
        .collect()
}

/// Validate that tetrahedral stereo has 3 or 4 neighbors.
pub(super) fn validate_tetrahedral_geometry(
    mol: &TableMolecule,
    atom_idx: usize,
) -> Result<(), RaiseError> {
    let count = neighbor_count(mol, atom_idx);
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
    bond_idx: usize,
) -> Result<(), RaiseError> {
    let bond = &mol.bonds[bond_idx];
    if bond.order != TableBondOrder::Single || bond.direction.is_none() {
        return Ok(());
    }
    let flanks_capable = [bond.start_atom() as usize, bond.end_atom() as usize]
        .into_iter()
        .any(|atom| {
            mol.bonds.iter().any(|d| {
                d.order == TableBondOrder::Double
                    && d.atoms
                        .other(atom as u32)
                        .is_some_and(|partner| cis_trans_capable(mol, atom, partner as usize))
            })
        });
    if flanks_capable {
        Ok(())
    } else {
        Err(RaiseError::DanglingBondDirection { bond: bond_idx })
    }
}

/// Double bond is cis-trans capable iff both ends have distinct substituents.
pub(super) fn cis_trans_capable(mol: &TableMolecule, atom_1: usize, atom_2: usize) -> bool {
    let side_1: Vec<_> = atom_ordering(mol, atom_1)
        .into_iter()
        .filter(|&atom| atom != atom_2)
        .collect();
    let side_2: Vec<_> = atom_ordering(mol, atom_2)
        .into_iter()
        .filter(|&atom| atom != atom_1)
        .collect();
    !side_1.is_empty() && !side_2.is_empty() && side_1.iter().all(|atom| !side_2.contains(atom))
}

/// Arrangement of the bond atom `atom_idx` of stereogenic double bond. Errors when its markers disagree.
pub(super) fn cis_trans_side(
    mol: &TableMolecule,
    atom_idx: usize,
    other_atom_idx: usize,
) -> Result<Option<StereoBondAtom>, RaiseError> {
    let substituents: Vec<usize> = atom_ordering(mol, atom_idx)
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
    let toward_first = direction(mol, atom_idx, first);
    let toward_second = substituents
        .get(1)
        .and_then(|&second| direction(mol, atom_idx, second))
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

/// Halfplane (top/bottom) of `other_atom_idx` viewed from `atom_idx`.
fn direction(
    mol: &TableMolecule,
    atom_idx: usize,
    other_atom_idx: usize,
) -> Option<StereoHalfplane> {
    mol.bonds.iter().find_map(|bond| {
        if bond.order != TableBondOrder::Single {
            return None;
        }
        if bond.atoms.other(atom_idx as u32)? as usize != other_atom_idx {
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
