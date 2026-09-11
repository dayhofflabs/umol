//! Wedge-based atom stereo interpretation and neighbor ordering for raise.

use umol_geometric_core::{complementary_direction, signed_volume, Point3D, AXIS_SIDE_TOLERANCE};
use umol_graph_ir::ir::NoncovalentBondKind;

use super::RaiseError;
use crate::table_ir::bond::{BondNoncovalent as TableNoncovalent, BondOrder as TableBondOrder};
use crate::table_ir::{AtomNeighbors, BondOrientation, Molecule as TableMolecule};

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

/// Out-of-plane direction.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum StereoOutofPlane {
    Front,
    Back,
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
