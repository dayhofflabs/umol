//! TableIR → `umol_ast::MoleculeAst` raise.
//!
//! Implements `TryIntoAst<MoleculeAst> for &Molecule` (and the per-atom and
//! per-bond analogues). Table IR fields copy to `Lit` / `Undetermined`; IO
//! raise applies fixed IO ground semantics for resolution.

use std::any::Any;
use std::collections::HashSet;

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraint, AtomId, BondAst, BondConstraint, Constraints,
    DativeBondAst, ElementAst, IsotopeMassAst, MoleculeAst, MulticenterBondAst, NoncovalentBondAst,
    NoncovalentBondKind, SpinStateAst, StereoConfigurationAst, StereoCosetAst, TryIntoAst,
    ValueAst,
};
use umol_geometric_core::{complementary_direction, signed_volume, Point3D};
use umol_perm::{space, ClassKey, Permutation};
use umol_shared::error::UmolError;

use crate::table_ir::atom::Atom as TableAtom;
use crate::table_ir::bond::{
    Bond as TableBond, BondDonation as TableBondDonation, BondNoncovalent as TableNoncovalent,
    BondOrder as TableBondOrder,
};
use crate::table_ir::{
    BondStereo, BondWedge, Chirality, ChiralityFrame, Molecule as TableMolecule,
};

/// Error variants for TableIR -> MoleculeAst raise.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RaiseError {
    /// A tetrahedral stereo assertion on an atom that cannot host a tetrahedral center.
    #[error("tetrahedral stereo at atom {atom} with {count} ligands, expect 3 or 4 ligands")]
    TetrahedralLigandCount { atom: usize, count: usize },
    /// A cis/trans double-bond atom with no substituent besides the other end.
    #[error(
        "cis/trans stereo at atom {atom} with no substituent besides the other double-bond end"
    )]
    CisTransLigandCount { atom: usize },
    /// Directional `/`,`\` bonds on one double-bond atom that imply both faces for its substituent.
    #[error("contradictory cis/trans markers at atom {atom}")]
    CisTransConflict { atom: usize },
    /// Multiple wedge bonds at a tetrahedral center that disagree on the configuration.
    #[error("inconsistent wedge bonds at atom {atom}")]
    WedgeConflict { atom: usize },
}

impl UmolError for RaiseError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// A ligand position in a raise-time stereo ordering: a neighbor atom, or an opaque virtual ligand
/// (an under-determined implicit H or lone pair — raise asserts the coset without deciding which).
/// Raise-local and distinct from `umol_ast::ast::StereoLigand`; `Permutation::between` is generic over
/// `Eq`, so the orderings carry no AST type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum StereoLigand {
    Atom(usize),
    /// An under-determined virtual ligand, tagged by its host atom so the two double-bond atoms'
    /// virtuals in a `#C` ordering stay distinct (`Permutation::between` matches ligands by equality).
    Virtual(usize),
}

impl TryIntoAst<MoleculeAst> for &TableMolecule {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ast(self, ctx: &Self::Ctx) -> Result<MoleculeAst, RaiseError> {
        let mut atoms: Vec<AtomAst> = self
            .atoms
            .iter()
            .map(|a| a.try_into_ast(ctx))
            .collect::<Result<_, _>>()?;

        let mut bonds = Vec::new();
        let mut dative_bonds: Vec<(Vec<AtomId>, AtomId, DativeBondAst)> = Vec::new();
        let mut noncovalent_bonds = Vec::new();
        for (bond_idx, b) in self.bonds.iter().enumerate() {
            let a_idx = AtomId(b.atoms.first());
            let b_idx = AtomId(b.atoms.second());
            if let Some(kind) = b.noncovalent.map(noncovalent_kind) {
                noncovalent_bonds.push((a_idx, b_idx, NoncovalentBondAst::from_kind(kind)));
            } else if let Some(donation) = b.donation {
                let (donor, acceptor) = match donation {
                    TableBondDonation::Donating => (a_idx, b_idx),
                    TableBondDonation::Accepting => (b_idx, a_idx),
                    _ => {
                        bonds.push((a_idx, b_idx, b.try_into_ast(ctx)?));
                        continue;
                    }
                };
                let dative_bond = DativeBondAst::new(raise_bond_order(b.order));
                dative_bonds.push((vec![donor], acceptor, dative_bond));
            } else {
                let mut bond_ast = b.try_into_ast(ctx)?;
                if let Some(constraint) = raise_cis_trans_stereo(self, bond_idx)? {
                    bond_ast.constraints.add(constraint);
                }
                bonds.push((a_idx, b_idx, bond_ast));
            }
        }

        let multicenter_bond: Vec<(Vec<AtomId>, MulticenterBondAst)> = self
            .multicenter_bonds
            .iter()
            .map(|mc| {
                let mut seen = HashSet::new();
                let atoms: Vec<AtomId> = mc
                    .all_atoms()
                    .into_iter()
                    .filter(|a| seen.insert(*a))
                    .map(AtomId)
                    .collect();
                let n = atoms.len();
                (
                    atoms,
                    MulticenterBondAst::new(vec![ValueAst::Undetermined; n]),
                )
            })
            .collect();

        for (atom_idx, (table_atom, atom)) in self.atoms.iter().zip(atoms.iter_mut()).enumerate() {
            match table_atom.aromatic {
                Some(true) => {
                    atom.constraints.add(AtomConstraint::AromaticValence(
                        AromaticValenceAst::Aromatic(ValueAst::Undetermined),
                    ));
                }
                Some(false) => {
                    atom.constraints.add(AtomConstraint::AromaticValence(
                        AromaticValenceAst::NotAromatic,
                    ));
                }
                None => {}
            }

            if let Some(constraint) = raise_tetrahedral_stereo(self, atom_idx)? {
                atom.constraints.add(constraint);
            }

            if matches!(atom.isotope_mass, IsotopeMassAst::Undetermined) {
                atom.isotope_mass = IsotopeMassAst::Natural;
            }
            if matches!(atom.charge, ValueAst::Undetermined) {
                atom.charge = ValueAst::Lit(0);
            }
            if matches!(atom.spin.unpaired, ValueAst::Undetermined) {
                atom.spin.unpaired = ValueAst::Lit(0);
            }
            atom.constraints.retain(|c| !c.is_undetermined());
        }
        let constraints = Constraints::new();

        Ok(MoleculeAst::from_parts(
            atoms,
            bonds,
            dative_bonds,
            vec![],
            multicenter_bond,
            noncovalent_bonds,
            Vec::new(),
            Vec::new(),
            constraints,
        ))
    }
}

impl TryIntoAst<AtomAst> for &TableAtom {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ast(self, _ctx: &Self::Ctx) -> Result<AtomAst, RaiseError> {
        Ok(AtomAst {
            element: ElementAst::Lit(self.element),
            isotope_mass: match self.isotope_mass {
                Some(m) => IsotopeMassAst::Lit(m as i64),
                None => IsotopeMassAst::Undetermined,
            },
            charge: match self.charge {
                Some(c) => ValueAst::Lit(c as i64),
                None => ValueAst::Undetermined,
            },
            implicit_hydrogens: match self.implicit_hydrogens {
                Some(n) => ValueAst::Lit(n as i64),
                None => ValueAst::Undetermined,
            },
            lone_pairs: match self.lone_pairs {
                Some(n) => ValueAst::Lit(n as i64),
                None => ValueAst::Undetermined,
            },
            spin: match self.unpaired_electrons {
                Some(u) => SpinStateAst {
                    unpaired: ValueAst::Lit(u as i64),
                    multiplicity: ValueAst::Undetermined,
                },
                None => SpinStateAst::default(),
            },
            constraints: Default::default(),
        })
    }
}

impl TryIntoAst<BondAst> for &TableBond {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ast(self, _ctx: &Self::Ctx) -> Result<BondAst, RaiseError> {
        let mut bond = BondAst::new(raise_bond_order(self.order));
        bond.charge = match self.charge {
            Some(c) => ValueAst::Lit(c as i64),
            None => ValueAst::Undetermined,
        };
        if matches!(self.order, TableBondOrder::Aromatic) {
            bond.constraints.add(BondConstraint::Aromatic);
        }
        Ok(bond)
    }
}

fn raise_bond_order(order: TableBondOrder) -> ValueAst {
    match order {
        TableBondOrder::Zero => ValueAst::Lit(0),
        TableBondOrder::Single => ValueAst::Lit(1),
        TableBondOrder::Double => ValueAst::Lit(2),
        TableBondOrder::Triple => ValueAst::Lit(3),
        TableBondOrder::Quadruple => ValueAst::Lit(4),
        TableBondOrder::Quintuple => ValueAst::Lit(5),
        TableBondOrder::Sextuple => ValueAst::Lit(6),
        // Definite-aromatic: localized bond order is 1 by Kekulé convention;
        // the aromatic flag is added separately as `BondConstraint::Aromatic`.
        // Renders as `1#a`.
        TableBondOrder::Aromatic => ValueAst::Lit(1),
        // Fuzzy orders: no concrete bond order can be assigned; raise to
        // `Undetermined`. Aromatic-flag setting (where applicable) is left
        // off — the chemistry of these is too ambiguous for the raise.
        // TODO: Convert to LitSet.
        TableBondOrder::SingleOrDouble
        | TableBondOrder::SingleOrAromatic
        | TableBondOrder::DoubleOrAromatic
        | TableBondOrder::Any => ValueAst::Undetermined,
    }
}

fn noncovalent_kind(kind: TableNoncovalent) -> NoncovalentBondKind {
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
fn tetrahedral_target_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
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
fn first_neighbor_toward_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
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
fn last_neighbor_away_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
    tetrahedral_target_ordering(mol, atom_idx)
}

/// The configuration index a wedge realizes in `ordering`: lift the wedged neighbor out of
/// plane (Up: +, Down: −), keep the rest in the 2D depiction (a virtual ligand at the
/// complementary in-plane direction), and take the sign of the signed volume of the four
/// points. The 0/1 mapping is pinned by the tests.
fn wedge_winding(
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
fn wedged_neighbors(mol: &TableMolecule, atom_idx: usize) -> Vec<(usize, bool)> {
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
fn tetrahedral_ligand_gate(mol: &TableMolecule, atom_idx: usize) -> Result<(), RaiseError> {
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

/// One double-bond atom's two #C ligands: by atom index (`first_ligand`, `second_ligand`; the second
/// is the under-determined `Virtual` when the atom has one substituent) and by face (`up_ligand`,
/// `down_ligand`).
struct StereoBondAtom {
    first_ligand: StereoLigand,
    second_ligand: StereoLigand,
    up_ligand: StereoLigand,
    down_ligand: StereoLigand,
}

/// The directional `/`,`\` of the single bond `atom_idx`–`other_atom_idx` as a face (`true` = up),
/// oriented toward `atom_idx` (inverted when `atom_idx` is the bond's end). `None` when the bond has
/// no `/`,`\`.
fn direction(mol: &TableMolecule, atom_idx: usize, other_atom_idx: usize) -> Option<bool> {
    mol.bonds.iter().find_map(|bond| {
        if bond.order != TableBondOrder::Single {
            return None;
        }
        if bond.atoms.other(atom_idx as u32)? as usize != other_atom_idx {
            return None;
        }
        let up = match bond.wedge? {
            BondWedge::Up => true,
            BondWedge::Down => false,
            _ => return None,
        };
        Some(if bond.start_atom() as usize == atom_idx {
            up
        } else {
            !up
        })
    })
}

/// One double-bond atom's #C side. `None` when it carries no `/`,`\` (coordinates / external).
/// Errors when it has no substituent besides `other_atom_idx`, or its markers disagree.
fn cis_trans_side(
    mol: &TableMolecule,
    atom_idx: usize,
    other_atom_idx: usize,
) -> Result<Option<StereoBondAtom>, RaiseError> {
    let substituents: Vec<usize> = atom_ordering(mol, atom_idx)
        .into_iter()
        .filter(|&n| n != other_atom_idx)
        .collect();
    let first = *substituents
        .first()
        .ok_or(RaiseError::CisTransLigandCount { atom: atom_idx })?;
    let first_ligand = StereoLigand::Atom(first);
    let second_ligand = substituents
        .get(1)
        .map_or(StereoLigand::Virtual(atom_idx), |&second| {
            StereoLigand::Atom(second)
        });
    // The first ligand's face: from the bond toward it and the inverted bond toward the second.
    let toward_first = direction(mol, atom_idx, first);
    let toward_second = substituents
        .get(1)
        .and_then(|&second| direction(mol, atom_idx, second))
        .map(|up| !up);
    let up = match (toward_first, toward_second) {
        (Some(a), Some(b)) if a != b => {
            return Err(RaiseError::CisTransConflict { atom: atom_idx })
        }
        (Some(face), _) | (_, Some(face)) => face,
        (None, None) => return Ok(None),
    };
    Ok(Some(if up {
        StereoBondAtom {
            first_ligand,
            second_ligand,
            up_ligand: first_ligand,
            down_ligand: second_ligand,
        }
    } else {
        StereoBondAtom {
            first_ligand,
            second_ligand,
            up_ligand: second_ligand,
            down_ligand: first_ligand,
        }
    }))
}

/// Raise tetrahedral stereo constraint for `atom_idx`. raise asserts only the `#T` coset — it does
/// not decide whether an under-determined fourth ligand is an implicit H or a lone pair, and writes
/// no atom fields. Errors when the atom cannot host a tetrahedral center; `Ok(None)` when there is
/// no tetrahedral descriptor.
fn raise_tetrahedral_stereo(
    mol: &TableMolecule,
    atom_idx: usize,
) -> Result<Option<AtomConstraint>, RaiseError> {
    let chirality = mol.atoms[atom_idx].chirality;
    let (source, target, source_idx): (Vec<StereoLigand>, Vec<StereoLigand>, usize) =
        match (chirality, mol.chirality_frame) {
            (Some(Chirality::Unspecified), _) => {
                tetrahedral_ligand_gate(mol, atom_idx)?;
                return Ok(Some(AtomConstraint::TetrahedralStereo(
                    StereoConfigurationAst::stereo(StereoCosetAst::Undetermined),
                )));
            }
            (Some(token), Some(ChiralityFrame::FirstNeighborToward)) => {
                let source_idx = match token {
                    Chirality::CounterClockwise | Chirality::Tetrahedral { arr: 1 } => 0,
                    Chirality::Clockwise | Chirality::Tetrahedral { arr: 2 } => 1,
                    _ => return Ok(None),
                };
                tetrahedral_ligand_gate(mol, atom_idx)?;
                (
                    first_neighbor_toward_ordering(mol, atom_idx),
                    tetrahedral_target_ordering(mol, atom_idx),
                    source_idx,
                )
            }
            (Some(token), Some(ChiralityFrame::LastNeighborAway)) => {
                let source_idx = match token {
                    Chirality::Clockwise => 0,
                    Chirality::CounterClockwise => 1,
                    _ => return Ok(None),
                };
                tetrahedral_ligand_gate(mol, atom_idx)?;
                (
                    last_neighbor_away_ordering(mol, atom_idx),
                    tetrahedral_target_ordering(mol, atom_idx),
                    source_idx,
                )
            }
            (None, _) => {
                let Some(positions) = mol.positions.as_ref() else {
                    return Ok(None);
                };
                let wedges = wedged_neighbors(mol, atom_idx);
                let Some(&(wedged, up)) = wedges.first() else {
                    return Ok(None);
                };
                tetrahedral_ligand_gate(mol, atom_idx)?;
                let ordering = tetrahedral_target_ordering(mol, atom_idx);
                let source_idx = wedge_winding(positions, atom_idx, &ordering, wedged, up);
                for &(other_wedged, other_up) in &wedges[1..] {
                    if wedge_winding(positions, atom_idx, &ordering, other_wedged, other_up)
                        != source_idx
                    {
                        return Err(RaiseError::WedgeConflict { atom: atom_idx });
                    }
                }
                (ordering.clone(), ordering, source_idx)
            }
            (Some(_), None) => return Ok(None),
        };
    let coset = space(ClassKey::Tetrahedral)
        .reindex(source_idx as u32, Permutation::between(&source, &target));
    Ok(Some(AtomConstraint::TetrahedralStereo(
        StereoConfigurationAst::stereo(StereoCosetAst::Lit(coset)),
    )))
}

/// Whether either double-bond atom carries a `/`,`\` on one of its substituent bonds — i.e. there
/// is a cis/trans assertion to interpret (otherwise the bond is plain or coordinate-only).
fn has_cis_trans_marker(mol: &TableMolecule, atom_1: usize, atom_2: usize) -> bool {
    [atom_1, atom_2].into_iter().any(|atom_idx| {
        atom_ordering(mol, atom_idx)
            .into_iter()
            .filter(|&other| other != atom_1 && other != atom_2)
            .any(|other| direction(mol, atom_idx, other).is_some())
    })
}

/// Raise cis/trans stereo constraint for `bond_idx` via the `#C` coset space (`D₄`/`V`): build each
/// side's by-face and atom-index pairs, then take the coset of the relabeling between them.
fn raise_cis_trans_stereo(
    mol: &TableMolecule,
    bond_idx: usize,
) -> Result<Option<BondConstraint>, RaiseError> {
    let bond = &mol.bonds[bond_idx];
    if bond.order != TableBondOrder::Double {
        return Ok(None);
    }
    if bond.stereo == Some(BondStereo::Either) {
        return Ok(Some(BondConstraint::CisTransStereo(
            StereoConfigurationAst::stereo(StereoCosetAst::Undetermined),
        )));
    }
    let atom_1 = bond.start_atom() as usize;
    let atom_2 = bond.end_atom() as usize;
    // No `/`,`\` anywhere on the flanking bonds: not a cis/trans assertion (plain or coordinates).
    if !has_cis_trans_marker(mol, atom_1, atom_2) {
        return Ok(None);
    }
    let (Some(side_1), Some(side_2)) = (
        cis_trans_side(mol, atom_1, atom_2)?,
        cis_trans_side(mol, atom_2, atom_1)?,
    ) else {
        return Ok(None);
    };
    let source = [
        side_1.up_ligand,
        side_1.down_ligand,
        side_2.up_ligand,
        side_2.down_ligand,
    ];
    let target = [
        side_1.first_ligand,
        side_1.second_ligand,
        side_2.first_ligand,
        side_2.second_ligand,
    ];
    let coset = space(ClassKey::CisTrans).index(Permutation::between(&source, &target));
    Ok(Some(BondConstraint::CisTransStereo(
        StereoConfigurationAst::stereo(StereoCosetAst::Lit(coset)),
    )))
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::BondId;
    use umol_shared::element::Element;

    use super::*;
    use crate::ctfile::parse_mol_to_ast;
    use crate::ctfile::parser::parse_mol_bytes_to_table_ir;
    use crate::smiles::parse_smiles_to_ast;
    use crate::smiles::parser::parse_smiles_bytes_to_table_ir;
    use crate::table_ir::atom::Atom as TableAtom;
    use crate::table_ir::bond::{Bond as TableBond, BondOrder as TableBondOrder};
    use crate::table_ir::Molecule as TableMolecule;

    #[fixture]
    fn methane() -> TableMolecule {
        let mut atom = TableAtom::from_element(Element::C);
        atom.implicit_hydrogens = Some(4);
        let mut mol = TableMolecule::empty();
        mol.atoms.push(atom);
        mol
    }

    #[fixture]
    fn carbon() -> TableMolecule {
        let mut mol = TableMolecule::empty();
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol
    }

    #[fixture]
    fn diatomic(#[default(TableBondOrder::Single)] order: TableBondOrder) -> TableMolecule {
        let mut mol = TableMolecule::empty();
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.bonds.push(TableBond::new(0, 1, order));
        mol
    }

    const METHANE_MOL: &str = "Methane\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const BENZENE_AROMATIC_MOL: &str = "benzene\n\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  4  0  0  0  0\n  2  3  4  0  0  0  0\n  3  4  4  0  0  0  0\n  4  5  4  0  0  0  0\n  5  6  4  0  0  0  0\n  6  1  4  0  0  0  0\nM  END\n";

    const CARBON_H0_EXPLICIT_MOL: &str = "carbon-h0\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  1  0  0  0  0  0  0  0  0\nM  END\n";

    const CHIRAL_PARITY_MOL: &str = "chiral\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n";

    const CIS_TRANS_EITHER_MOL: &str = "butene\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  3  0  0  0\n  3  4  1  0  0  0  0\nM  END\n";

    #[rstest]
    fn test_table_molecule_try_into_ast(methane: TableMolecule) {
        let ast: MoleculeAst = (&methane).try_into_ast(&()).unwrap();
        assert_eq!(ast.atoms().count(), 1);
        let atom = ast.atom(AtomId(0)).ast;
        assert_eq!(atom.element, ElementAst::Lit(Element::C));
        assert!(matches!(atom.implicit_hydrogens, ValueAst::Lit(4)));
        assert!(matches!(atom.charge, ValueAst::Lit(0)));
        assert!(matches!(atom.lone_pairs, ValueAst::Undetermined));
        assert!(matches!(atom.spin.unpaired, ValueAst::Lit(0)));
    }

    #[rstest]
    #[case::table_aromatic_none(None, AromaticValenceAst::Undetermined)]
    #[case::table_aromatic_false(Some(false), AromaticValenceAst::NotAromatic)]
    #[case::table_aromatic_true(Some(true), AromaticValenceAst::Aromatic(ValueAst::Undetermined))]
    fn test_table_molecule_try_into_ast_aromatic(
        mut carbon: TableMolecule,
        #[case] aromatic: Option<bool>,
        #[case] expected: AromaticValenceAst,
    ) {
        carbon.atoms[0].aromatic = aromatic;
        let ast: MoleculeAst = (&carbon).try_into_ast(&()).unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.constraints.aromatic_valence(),
            expected
        );
    }

    #[rstest]
    fn test_table_molecule_try_into_ast_bond_order(
        #[with(TableBondOrder::Double)] diatomic: TableMolecule,
    ) {
        let ast: MoleculeAst = (&diatomic).try_into_ast(&()).unwrap();
        let bond = ast.bond(BondId(0)).ast;
        assert!(matches!(bond.order, ValueAst::Lit(2)));
    }

    #[rstest]
    fn test_table_molecule_try_into_ast_aromatic_bond(
        #[with(TableBondOrder::Aromatic)] diatomic: TableMolecule,
    ) {
        let ast: MoleculeAst = (&diatomic).try_into_ast(&()).unwrap();
        let bond = ast.bond(BondId(0)).ast;
        assert!(matches!(bond.order, ValueAst::Lit(1)));
        assert!(bond
            .constraints
            .iter()
            .any(|c| matches!(c, BondConstraint::Aromatic)));
        for i in 0..2 {
            assert!(matches!(
                ast.atom(AtomId(i)).ast.constraints.aromatic_valence(),
                AromaticValenceAst::Undetermined
            ));
        }
    }

    #[rstest]
    #[case::methane(METHANE_MOL, "C#i=#c0#u0")]
    #[case::benzene(BENZENE_AROMATIC_MOL, "C#i=#c0#u0")]
    #[case::carbon_h0(CARBON_H0_EXPLICIT_MOL, "C#i=#c0#h0#u0")]
    fn test_parse_mol_to_ast(#[case] input: &str, #[case] expected_atom: &str) {
        let ast = parse_mol_to_ast(input).unwrap();
        let atom = ast.atom(AtomId(0)).ast;
        assert_eq!(atom.charge, ValueAst::Lit(0));
        assert!(matches!(
            atom.constraints.aromatic_valence(),
            AromaticValenceAst::Undetermined
        ));
        assert_eq!(atom.to_string(), expected_atom);
    }

    #[rstest]
    #[case::organic("C", "C#i=#c0#u0#a!")]
    fn test_parse_smiles_to_ast(#[case] input: &str, #[case] expected_atom: &str) {
        let ast = parse_smiles_to_ast(input).unwrap();
        let atom = ast.atom(AtomId(0)).ast;
        assert_eq!(atom.charge, ValueAst::Lit(0));
        assert!(matches!(atom.implicit_hydrogens, ValueAst::Undetermined));
        assert!(matches!(
            atom.constraints.aromatic_valence(),
            AromaticValenceAst::NotAromatic
        ));
        assert_eq!(atom.to_string(), expected_atom);
    }

    #[rstest]
    #[case::cfclbri_clockwise(parse_smiles_bytes_to_table_ir(b"Br[C@@](F)(Cl)I").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::cfclbri_counterclockwise(parse_smiles_bytes_to_table_ir(b"Br[C@](F)(Cl)I").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::cfclbri_fluorine_first(parse_smiles_bytes_to_table_ir(b"F[C@](Cl)(Br)I").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::methyloxirane_explicit_h(parse_smiles_bytes_to_table_ir(b"C[C@@]1([H])OC1").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::butan_2_ol(parse_smiles_bytes_to_table_ir(b"C[C@@H](O)CC").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::alanine(parse_smiles_bytes_to_table_ir(b"C[C@H](N)C(O)=O").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::ring_then_branch(parse_smiles_bytes_to_table_ir(b"C[C@]1(Cl)CC(C)CC1").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::branch_then_ring(parse_smiles_bytes_to_table_ir(b"C[C@](Cl)1CC(C)CC1").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::mol_parity_clockwise(parse_mol_bytes_to_table_ir(CHIRAL_PARITY_MOL.as_bytes()).unwrap(), 0, Some(StereoCosetAst::Lit(0)))]
    #[case::sulfoxide_counterclockwise(parse_smiles_bytes_to_table_ir(b"C[S@](=O)CC").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::sulfoxide_clockwise(parse_smiles_bytes_to_table_ir(b"C[S@@](=O)CC").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::sulfoxide_charge_separated(parse_smiles_bytes_to_table_ir(b"C[S@@+]([O-])CC").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::no_descriptor(parse_smiles_bytes_to_table_ir(b"F[C@](Cl)(Br)I").unwrap(), 0, None)]
    fn test_raise_tetrahedral_stereo(
        #[case] mol: TableMolecule,
        #[case] atom_idx: usize,
        #[case] expected: Option<StereoCosetAst>,
    ) {
        let expected = expected
            .map(|coset| AtomConstraint::TetrahedralStereo(StereoConfigurationAst::stereo(coset)));
        assert_eq!(raise_tetrahedral_stereo(&mol, atom_idx), Ok(expected));
    }

    #[rstest]
    #[case::dimethyl_sulfide(parse_smiles_bytes_to_table_ir(b"C[S@]C").unwrap(), 1, 2)]
    fn test_raise_tetrahedral_stereo_error(
        #[case] mol: TableMolecule,
        #[case] atom_idx: usize,
        #[case] count: usize,
    ) {
        assert_eq!(
            raise_tetrahedral_stereo(&mol, atom_idx),
            Err(RaiseError::TetrahedralLigandCount {
                atom: atom_idx,
                count,
            })
        );
    }

    #[rstest]
    #[case::trans(parse_smiles_bytes_to_table_ir(b"F/C=C/F").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::cis(parse_smiles_bytes_to_table_ir(b"F/C=C\\F").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::fluoropropene_e(parse_smiles_bytes_to_table_ir(b"F/C=C/C").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::fluoropropene_z(parse_smiles_bytes_to_table_ir(b"F/C=C\\C").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::fluoropropene_z_flipped(parse_smiles_bytes_to_table_ir(b"F\\C=C/C").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::fluoropropene_z_methyl_first(parse_smiles_bytes_to_table_ir(b"C/C=C\\F").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::azomethane_e(parse_smiles_bytes_to_table_ir(b"C/N=N/C").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::azomethane_z(parse_smiles_bytes_to_table_ir(b"C/N=N\\C").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::ethylideneoxirane(parse_smiles_bytes_to_table_ir(b"C/C=C1CO\\1").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::butanone_oxime(parse_smiles_bytes_to_table_ir(b"C/C(CC)=N\\O").unwrap(), 3, Some(StereoCosetAst::Lit(0)))]
    #[case::fluoropropene_e_backslash(parse_smiles_bytes_to_table_ir(b"F\\C=C\\C").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::fluoropropene_e_methyl_first(parse_smiles_bytes_to_table_ir(b"C/C=C/F").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::fluoropropene_e_methyl_first_backslash(parse_smiles_bytes_to_table_ir(b"C\\C=C\\F").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::trisubstituted(parse_smiles_bytes_to_table_ir(b"F/C(C)=C(Cl)/C").unwrap(), 2, Some(StereoCosetAst::Lit(0)))]
    #[case::mol_either(parse_mol_bytes_to_table_ir(CIS_TRANS_EITHER_MOL.as_bytes()).unwrap(), 1, Some(StereoCosetAst::Undetermined))]
    #[case::one_sided_marker(parse_smiles_bytes_to_table_ir(b"C(C)=C(Cl)/C").unwrap(), 1, None)]
    #[case::plain_double(parse_smiles_bytes_to_table_ir(b"C=C").unwrap(), 0, None)]
    fn test_raise_cis_trans_stereo(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: Option<StereoCosetAst>,
    ) {
        let expected = expected
            .map(|coset| BondConstraint::CisTransStereo(StereoConfigurationAst::stereo(coset)));
        assert_eq!(raise_cis_trans_stereo(&mol, bond_idx), Ok(expected));
    }

    #[rstest]
    #[case::conflict(parse_smiles_bytes_to_table_ir(b"F/C(\\Cl)=CF").unwrap(), 2, RaiseError::CisTransConflict { atom: 1 })]
    #[case::no_substituent(parse_smiles_bytes_to_table_ir(b"F/C=C").unwrap(), 1, RaiseError::CisTransLigandCount { atom: 2 })]
    fn test_raise_cis_trans_stereo_error(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: RaiseError,
    ) {
        assert_eq!(raise_cis_trans_stereo(&mol, bond_idx), Err(expected));
    }
}
