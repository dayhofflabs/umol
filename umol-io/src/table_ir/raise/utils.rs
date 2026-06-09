//! Raise-time stereo helpers for `super`: ligand orderings, directional-bond faces, and the
//! capability/validation predicates the `raise_*` functions build on.

use umol_ast::ast::NoncovalentBondKind;
use umol_geometric_core::{complementary_direction, signed_volume, Point3D};

use crate::table_ir::bond::{BondNoncovalent as TableNoncovalent, BondOrder as TableBondOrder};
use crate::table_ir::{BondWedge, Molecule as TableMolecule};

use super::RaiseError;

/// A ligand position in a raise-time stereo ordering: a neighbor atom, or an opaque virtual ligand
/// (an under-determined implicit H or lone pair — raise asserts the coset without deciding which).
/// Raise-local and distinct from `umol_ast::ast::StereoLigand`; `Permutation::between` is generic over
/// `Eq`, so the orderings carry no AST type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum StereoLigand {
    Atom(usize),
    /// An under-determined virtual ligand, tagged by its host atom so the two double-bond atoms'
    /// virtuals in a `#C` ordering stay distinct (`Permutation::between` matches ligands by equality).
    Virtual(usize),
}

pub(super) fn noncovalent_kind(kind: TableNoncovalent) -> NoncovalentBondKind {
    match kind {
        TableNoncovalent::Hydrogen => NoncovalentBondKind::HydrogenBond,
    }
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
    indices
}

/// Neighbor atom ordering of `atom_idx` by bond ordering (used by SMILES, which refers to it as parse ordering).
/// Neighbor atoms appear in the order of incident bonds (including ring-closure indices).
/// Can differ from the ascending atom index order when rings are present.
fn bond_neighbor_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<usize> {
    mol.bonds
        .iter()
        .filter_map(|bond| bond.atoms.other(atom_idx as u32))
        .map(|other| other as usize)
        .collect()
}

/// #T umol target ordering: neighbors ascending as `Atom`, then one `Virtual` (the under-determined
/// fourth ligand, placed last) when the center has three neighbors.
pub(super) fn tetrahedral_target_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
    let mut ordering: Vec<StereoLigand> = atom_ordering(mol, atom_idx)
        .into_iter()
        .map(StereoLigand::Atom)
        .collect();
    if ordering.len() == 3 {
        ordering.push(StereoLigand::Virtual(atom_idx));
    }
    ordering
}

/// #T source ordering, FirstNeighborToward (SMILES/SMARTS): neighbors in parse order, with the
/// under-determined `Virtual` at the implicit slot (first if `atom_idx` opens the SMILES, else second)
/// when the center has three neighbors.
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

/// #T source ordering, LastNeighborAway (MDL / MOL parity). Atom-number order = atom-index order with
/// the under-determined ligand last — identical to the target frame, so it passes through unchanged;
/// only the parity winding (read with that last ligand behind) converts.
pub(super) fn last_neighbor_away_ordering(
    mol: &TableMolecule,
    atom_idx: usize,
) -> Vec<StereoLigand> {
    tetrahedral_target_ordering(mol, atom_idx)
}

/// The configuration index a wedge realizes in `ordering`: lift the wedged neighbor out of
/// plane (Up: +, Down: −), keep the rest in the 2D depiction (a virtual ligand at the
/// complementary in-plane direction), and take the sign of the signed volume of the four
/// points. The 0/1 mapping is pinned by the tests.
pub(super) fn wedge_winding(
    positions: &[Point3D],
    atom_idx: usize,
    ordering: &[StereoLigand],
    wedged: usize,
    up: bool,
) -> usize {
    let z = if up { 1.0 } else { -1.0 };
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
    if signed_volume(points[0], points[1], points[2], points[3]) > 0.0 {
        0
    } else {
        1
    }
}

/// The neighbors reached by wedge bonds at `atom_idx`, each with its up/down sense (`true` = up).
pub(super) fn wedged_neighbors(mol: &TableMolecule, atom_idx: usize) -> Vec<(usize, bool)> {
    mol.bonds
        .iter()
        .filter_map(
            |bond| match (bond.atoms.other(atom_idx as u32), bond.wedge) {
                (Some(other), Some(BondWedge::Up)) => Some((other as usize, true)),
                (Some(other), Some(BondWedge::Down)) => Some((other as usize, false)),
                _ => None,
            },
        )
        .collect()
}

/// A tetrahedral center has four ligands, so it needs three neighbors (one under-determined virtual
/// completes it) or four. Fewer (or more) cannot host a tetrahedral assertion.
pub(super) fn tetrahedral_ligand_gate(mol: &TableMolecule, atom_idx: usize) -> Result<(), RaiseError> {
    let count = atom_ordering(mol, atom_idx).len();
    if count == 3 || count == 4 {
        Ok(())
    } else {
        Err(RaiseError::TetrahedralLigandCount {
            atom: atom_idx,
            count,
        })
    }
}

/// Which face of the double-bond axis a directional `/`,`\` puts a substituent on, as seen from a
/// given double-bond atom.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum StereoFace {
    Above,
    Below,
}

impl StereoFace {
    fn flip(self) -> Self {
        match self {
            Self::Above => Self::Below,
            Self::Below => Self::Above,
        }
    }
}

/// One double-bond atom's two #C ligands by atom index (`first_ligand`, `second_ligand`; the second is
/// the under-determined `Virtual` when the atom has one substituent), with the face of the first. The
/// second ligand sits on `first_face.flip()`.
pub(super) struct StereoBondAtom {
    pub(super) first_ligand: StereoLigand,
    pub(super) second_ligand: StereoLigand,
    pub(super) first_face: StereoFace,
}

/// The face the directional `/`,`\` of the single bond `atom_idx`–`other_atom_idx` puts
/// `other_atom_idx` on, seen from `atom_idx` (flipped when `atom_idx` is the bond's stored end).
/// `None` when the bond has no `/`,`\`.
fn direction(mol: &TableMolecule, atom_idx: usize, other_atom_idx: usize) -> Option<StereoFace> {
    mol.bonds.iter().find_map(|bond| {
        if bond.order != TableBondOrder::Single {
            return None;
        }
        if bond.atoms.other(atom_idx as u32)? as usize != other_atom_idx {
            return None;
        }
        let face = match bond.wedge? {
            BondWedge::Up => StereoFace::Above,
            BondWedge::Down => StereoFace::Below,
            _ => return None,
        };
        Some(if bond.start_atom() as usize == atom_idx {
            face
        } else {
            face.flip()
        })
    })
}

/// One double-bond atom's #C side. `None` when it has no substituent besides `other_atom_idx` (caller
/// gates this with `cis_trans_capable`) or carries no `/`,`\`. Errors when its markers disagree.
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
        .map(StereoFace::flip);
    let first_face = match (toward_first, toward_second) {
        (Some(a), Some(b)) if a != b => {
            return Err(RaiseError::CisTransConflict { atom: atom_idx })
        }
        (Some(face), _) | (_, Some(face)) => face,
        (None, None) => return Ok(None),
    };
    Ok(Some(StereoBondAtom {
        first_ligand,
        second_ligand,
        first_face,
    }))
}

/// A double bond `atom_1`=`atom_2` is cis/trans-capable when both ends carry a substituent besides
/// each other.
pub(super) fn cis_trans_capable(mol: &TableMolecule, atom_1: usize, atom_2: usize) -> bool {
    atom_ordering(mol, atom_1).into_iter().any(|n| n != atom_2)
        && atom_ordering(mol, atom_2).into_iter().any(|n| n != atom_1)
}

/// A directional `/`,`\` single bond must flank a cis/trans-capable double bond — one of its atoms is
/// in a double bond whose both ends are substituted. A marker flanking none is dangling, the cis/trans
/// analog of a chirality token on an under-coordinated atom. `Ok(())` for any non-directional bond.
pub(super) fn validate_bond_wedge(mol: &TableMolecule, bond_idx: usize) -> Result<(), RaiseError> {
    let bond = &mol.bonds[bond_idx];
    if bond.order != TableBondOrder::Single
        || !matches!(bond.wedge, Some(BondWedge::Up | BondWedge::Down))
    {
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
        Err(RaiseError::DanglingBondWedge { bond: bond_idx })
    }
}
