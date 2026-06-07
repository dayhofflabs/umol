//! TableIR → `umol_ast::MoleculeAst` raise.
//!
//! Implements `TryIntoAst<MoleculeAst> for &Molecule` (and the per-atom and
//! per-bond analogues). Table IR fields copy to `Lit` / `Undetermined`; IO
//! raise applies fixed IO ground semantics for resolution.

use std::collections::HashSet;

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraint, AtomId, BondAst, BondConstraint, Constraints,
    DativeBondAst, ElementAst, IsotopeMassAst, MoleculeAst, MulticenterBondAst, NoncovalentBondAst,
    NoncovalentBondKind, SpinStateAst, StereoConfigurationAst, StereoCosetAst, StereoLigand,
    StereoLigandKind, TryIntoAst, ValueAst,
};
use umol_perm::{space, ClassKey, Permutation};

use crate::position::Point3D;
use crate::table_ir::atom::Atom as TableAtom;
use crate::table_ir::bond::{
    Bond as TableBond, BondDonation as TableBondDonation, BondNoncovalent as TableNoncovalent,
    BondOrder as TableBondOrder,
};
use crate::table_ir::{BondStereo, BondWedge, Chirality, ChiralityFrame};
use crate::table_ir::Molecule as TableMolecule;

/// Error variants for TableIR -> MoleculeAst raise.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RaiseError {}

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
                if let Some(constraint) = raise_cis_trans_stereo(self, bond_idx) {
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

        for (atom_idx, (table_atom, atom)) in
            self.atoms.iter().zip(atoms.iter_mut()).enumerate()
        {
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

            if let Some(constraint) = raise_tetrahedral_stereo(self, atom_idx) {
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

/// Neighbors of `atom_idx`, ascending index (= atom-list / MDL atom-number order).
fn neighbors(mol: &TableMolecule, atom_idx: usize) -> Vec<usize> {
    let mut indices: Vec<usize> = mol
        .bonds
        .iter()
        .filter_map(|bond| bond.atoms.other(atom_idx as u32))
        .map(|other| other as usize)
        .collect();
    indices.sort_unstable();
    indices
}

/// `atom_idx`'s virtual ligands: implicit H then lone pairs, each hosted by `atom_idx`. Shared by the
/// target and source orderings so they hold the identical set (`Permutation::between` requires it).
fn virtual_ligands(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
    let atom = &mol.atoms[atom_idx];
    let host = AtomId(atom_idx as u32);
    let hydrogens = atom.implicit_hydrogens.unwrap_or(0);
    let lone_pairs = atom.lone_pairs.unwrap_or(0);
    (0..hydrogens)
        .map(|_| StereoLigand::new(host, StereoLigandKind::ImplicitHydrogen))
        .chain((0..lone_pairs).map(|_| StereoLigand::new(host, StereoLigandKind::LonePair)))
        .collect()
}

/// Neighbors of `atom_idx` in input (parse) order — the OpenSMILES write order, read off the
/// incident bonds in `mol.bonds` order. (Ring-closure slotting is pinned by the raise tests.)
fn input_neighbor_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<usize> {
    mol.bonds
        .iter()
        .filter_map(|bond| bond.atoms.other(atom_idx as u32))
        .map(|other| other as usize)
        .collect()
}

/// #T umol target ordering: real neighbors ascending, then the virtual ligands.
fn tetrahedral_target_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
    let mut ordering: Vec<StereoLigand> = neighbors(mol, atom_idx)
        .iter()
        .map(|&n| StereoLigand::new(AtomId(n as u32), StereoLigandKind::Atom))
        .collect();
    ordering.extend(virtual_ligands(mol, atom_idx));
    ordering
}

/// #T source ordering, FirstNeighborToward (SMILES/SMARTS): neighbors in input order, the implicit H
/// at slot 0 if `atom_idx` opened the SMILES (`atom_idx == 0`) else slot 1; remaining virtuals after.
fn first_neighbor_toward_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
    let mut ordering: Vec<StereoLigand> = input_neighbor_ordering(mol, atom_idx)
        .iter()
        .map(|&n| StereoLigand::new(AtomId(n as u32), StereoLigandKind::Atom))
        .collect();
    let mut virtuals = virtual_ligands(mol, atom_idx).into_iter();
    if mol.atoms[atom_idx].implicit_hydrogens.unwrap_or(0) > 0 {
        let slot = if atom_idx > 0 { 1 } else { 0 };
        ordering.insert(slot, virtuals.next().expect("implicit hydrogen present"));
    }
    ordering.extend(virtuals);
    ordering
}

/// #T source ordering, LastNeighborAway (MDL / MOL parity). Atom-number order = atom-index order with
/// the H/lone pair last — identical to the target frame, so it passes through unchanged; only the
/// parity winding (read with that last ligand behind) converts.
fn last_neighbor_away_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
    tetrahedral_target_ordering(mol, atom_idx)
}

/// The empty in-plane direction at `atom_idx`: opposite the sum of unit vectors to its real
/// neighbors — where a virtual ligand sits when a wedge is read against the 2D depiction.
fn leftover_inplane(positions: &[Point3D], atom_idx: usize, neighbor_indices: &[usize]) -> [f64; 3] {
    let center = positions[atom_idx];
    let (mut sum_x, mut sum_y) = (0.0, 0.0);
    for &neighbor in neighbor_indices {
        let (dx, dy) = (positions[neighbor].x - center.x, positions[neighbor].y - center.y);
        let length = (dx * dx + dy * dy).sqrt();
        if length > 0.0 {
            sum_x += dx / length;
            sum_y += dy / length;
        }
    }
    [center.x - sum_x, center.y - sum_y, 0.0]
}

/// The configuration index a wedge realizes in `ordering`: lift the wedged neighbor out of plane
/// (Up: +, Down: −), keep the rest in the 2D depiction (a virtual ligand at the leftover direction),
/// and take the sign of the signed winding (3×3 determinant). The 0/1 mapping is pinned by the tests.
fn wedge_winding(
    positions: &[Point3D],
    atom_idx: usize,
    ordering: &[StereoLigand],
    wedged: usize,
    up: bool,
) -> usize {
    let lift = if up { 1.0 } else { -1.0 };
    let real: Vec<usize> = ordering
        .iter()
        .filter(|ligand| ligand.kind == StereoLigandKind::Atom)
        .map(|ligand| ligand.atom_id.0 as usize)
        .collect();
    let leftover = leftover_inplane(positions, atom_idx, &real);
    let points: Vec<[f64; 3]> = ordering
        .iter()
        .map(|ligand| {
            let index = ligand.atom_id.0 as usize;
            match ligand.kind {
                StereoLigandKind::Atom if index == wedged => {
                    [positions[wedged].x, positions[wedged].y, lift]
                }
                StereoLigandKind::Atom => [positions[index].x, positions[index].y, 0.0],
                _ => leftover,
            }
        })
        .collect();
    let edge = |i: usize| {
        [
            points[i][0] - points[0][0],
            points[i][1] - points[0][1],
            points[i][2] - points[0][2],
        ]
    };
    let (u, v, w) = (edge(1), edge(2), edge(3));
    let determinant = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0]);
    if determinant > 0.0 {
        0
    } else {
        1
    }
}

/// #T source ordering from the wedge + 2D depiction (the wedge is consumed). Frame-free, so the
/// ordering returned is the target frame and the index carries the configuration. `None` if no wedge.
fn tetrahedral_wedge_ordering(
    mol: &TableMolecule,
    atom_idx: usize,
) -> Option<(Vec<StereoLigand>, usize)> {
    let positions = mol.positions.as_ref()?;
    let (wedged, up) =
        mol.bonds
            .iter()
            .find_map(|bond| match (bond.atoms.other(atom_idx as u32), bond.wedge) {
                (Some(other), Some(BondWedge::Up)) => Some((other as usize, true)),
                (Some(other), Some(BondWedge::Down)) => Some((other as usize, false)),
                _ => None,
            })?;
    let ordering = tetrahedral_target_ordering(mol, atom_idx);
    let source_idx = wedge_winding(positions, atom_idx, &ordering, wedged, up);
    Some((ordering, source_idx))
}

/// One sp² carbon's substituents (ascending index, then its virtual ligands).
fn cis_trans_side_ligands(mol: &TableMolecule, carbon: usize, partner: usize) -> Vec<StereoLigand> {
    let mut ordering: Vec<StereoLigand> = neighbors(mol, carbon)
        .into_iter()
        .filter(|&n| n != partner)
        .map(|n| StereoLigand::new(AtomId(n as u32), StereoLigandKind::Atom))
        .collect();
    ordering.extend(virtual_ligands(mol, carbon));
    ordering
}

/// #C umol target ordering: atom_1's side then atom_2's side.
fn cis_trans_target_ordering(
    mol: &TableMolecule,
    atom_1: usize,
    atom_2: usize,
) -> Vec<StereoLigand> {
    let mut ordering = cis_trans_side_ligands(mol, atom_1, atom_2);
    ordering.extend(cis_trans_side_ligands(mol, atom_2, atom_1));
    ordering
}

/// #C source ordering from the directional `/`,`\` flanking bonds. Frame-free, so the ordering
/// returned is the target frame and the index carries cis/trans (syn ⇒ 0, anti ⇒ 1). `None` when a
/// side has no directional bond (⇒ MOL field 0 ⇒ coordinates ⇒ external).
fn cis_trans_wedge_ordering(
    mol: &TableMolecule,
    atom_1: usize,
    atom_2: usize,
) -> Option<(Vec<StereoLigand>, usize)> {
    let side = |carbon: usize, partner: usize| -> Option<i8> {
        mol.bonds.iter().find_map(|single| {
            if single.order != TableBondOrder::Single {
                return None;
            }
            let other = single.atoms.other(carbon as u32)? as usize;
            if other == partner {
                return None;
            }
            let direction = match single.wedge? {
                BondWedge::Up => 1,
                BondWedge::Down => -1,
                _ => return None,
            };
            Some(if single.start_atom() as usize == carbon {
                direction
            } else {
                -direction
            })
        })
    };
    let direction_1 = side(atom_1, atom_2)?;
    let direction_2 = side(atom_2, atom_1)?;
    let source_idx = usize::from(direction_1 != direction_2);
    Some((cis_trans_target_ordering(mol, atom_1, atom_2), source_idx))
}

/// B4 — tetrahedral atom → `#T`.
fn raise_tetrahedral_stereo(mol: &TableMolecule, atom_idx: usize) -> Option<AtomConstraint> {
    let (source, source_idx): (Vec<StereoLigand>, usize) =
        match (mol.atoms[atom_idx].chirality, mol.chirality_frame) {
            (Some(Chirality::Unspecified), _) => {
                return Some(AtomConstraint::TetrahedralStereo(
                    StereoConfigurationAst::stereo(StereoCosetAst::Undetermined),
                ));
            }
            (Some(token), Some(ChiralityFrame::FirstNeighborToward)) => {
                let source_idx = match token {
                    Chirality::CounterClockwise | Chirality::Tetrahedral { arr: 1 } => 0,
                    Chirality::Clockwise | Chirality::Tetrahedral { arr: 2 } => 1,
                    _ => return None,
                };
                (first_neighbor_toward_ordering(mol, atom_idx), source_idx)
            }
            (Some(token), Some(ChiralityFrame::LastNeighborAway)) => {
                let source_idx = match token {
                    Chirality::Clockwise => 0,
                    Chirality::CounterClockwise => 1,
                    _ => return None,
                };
                (last_neighbor_away_ordering(mol, atom_idx), source_idx)
            }
            (None, _) => tetrahedral_wedge_ordering(mol, atom_idx)?,
            (Some(_), None) => return None,
        };
    let target = tetrahedral_target_ordering(mol, atom_idx);
    let coset = space(ClassKey::Tetrahedral)
        .reindex(source_idx as u32, Permutation::between(&source, &target));
    Some(AtomConstraint::TetrahedralStereo(
        StereoConfigurationAst::stereo(StereoCosetAst::Lit(coset)),
    ))
}

/// B5 — double bond → `#C`.
fn raise_cis_trans_stereo(mol: &TableMolecule, bond_idx: usize) -> Option<BondConstraint> {
    let bond = &mol.bonds[bond_idx];
    if bond.order != TableBondOrder::Double {
        return None;
    }
    if bond.stereo == Some(BondStereo::Either) {
        return Some(BondConstraint::CisTransStereo(
            StereoConfigurationAst::stereo(StereoCosetAst::Undetermined),
        ));
    }
    let atom_1 = bond.start_atom() as usize;
    let atom_2 = bond.end_atom() as usize;
    let (source, source_idx) = cis_trans_wedge_ordering(mol, atom_1, atom_2)?;
    let target = cis_trans_target_ordering(mol, atom_1, atom_2);
    let coset = space(ClassKey::CisTrans)
        .reindex(source_idx as u32, Permutation::between(&source, &target));
    Some(BondConstraint::CisTransStereo(
        StereoConfigurationAst::stereo(StereoCosetAst::Lit(coset)),
    ))
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::*;
    use umol_ast::ast::BondId;
    use umol_shared::element::Element;

    use super::*;
    use crate::io::ctfile::config::CtfileIoConfig;
    use crate::io::ctfile::{parse_mol_bytes_with, parse_mol_to_ast};
    use crate::io::smiles::parse_smiles_to_ast;
    use crate::io::smiles::parser::parse_smiles_bytes_to_table_ir;
    use crate::ops::model::{
        AromaticityModel, ChemistryModel, CountsModel, ElementScope, RingLimits, ValenceModel,
    };
    use crate::ops::valence::{CountsValence, ValenceTable};
    use crate::table_ir::atom::Atom as TableAtom;
    use crate::table_ir::bond::{Bond as TableBond, BondOrder as TableBondOrder};
    use crate::table_ir::Molecule as TableMolecule;

    fn methane_table() -> TableMolecule {
        let mut atom = TableAtom::from_element(Element::C);
        atom.implicit_hydrogens = Some(4);
        let mut mol = TableMolecule::empty();
        mol.atoms.push(atom);
        mol
    }

    const METHANE_MOL: &str = "Methane\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const BENZENE_AROMATIC_MOL: &str = "benzene\n\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  4  0  0  0  0\n  2  3  4  0  0  0  0\n  3  4  4  0  0  0  0\n  4  5  4  0  0  0  0\n  5  6  4  0  0  0  0\n  6  1  4  0  0  0  0\nM  END\n";

    const CARBON_H0_EXPLICIT_MOL: &str = "carbon-h0\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  1  0  0  0  0  0  0  0  0\nM  END\n";

    fn counts_model() -> CountsModel {
        CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        }
    }

    #[rstest]
    fn test_table_molecule_try_into_ast_methane_table() {
        let mol = methane_table();
        let ast: MoleculeAst = (&mol).try_into_ast(&()).unwrap();
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
    fn test_table_molecule_try_into_ast(
        #[case] table_aromatic: Option<bool>,
        #[case] expected: AromaticValenceAst,
    ) {
        let mut mol = TableMolecule::empty();
        let mut atom = TableAtom::from_element(Element::C);
        atom.aromatic = table_aromatic;
        mol.atoms.push(atom);
        let ast: MoleculeAst = (&mol).try_into_ast(&()).unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.constraints.aromatic_valence(),
            expected
        );
    }

    #[rstest]
    fn test_table_molecule_try_into_ast_bond_order() {
        let mut mol = TableMolecule::empty();
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.bonds.push(TableBond::new(0, 1, TableBondOrder::Double));
        let ast: MoleculeAst = (&mol).try_into_ast(&()).unwrap();
        let bond = ast.bond(BondId(0)).ast;
        assert!(matches!(bond.order, ValueAst::Lit(2)));
    }

    #[rstest]
    fn test_table_molecule_try_into_ast_aromatic_bond() {
        let mut mol = TableMolecule::empty();
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.atoms.push(TableAtom::from_element(Element::C));
        mol.bonds
            .push(TableBond::new(0, 1, TableBondOrder::Aromatic));
        let ast: MoleculeAst = (&mol).try_into_ast(&()).unwrap();
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
    #[case::methane(METHANE_MOL, 1, "C#i=#c0#h4#n0#u0#s#v0#a!")]
    #[case::benzene(BENZENE_AROMATIC_MOL, 6, "C#i=#c0#h#n0#u0#s#v2#a")]
    fn test_parse_mol_to_ast_counts_resolve(
        #[case] input: &str,
        #[case] atom_count: u32,
        #[case] expected_atom: &str,
    ) {
        let mut ast = parse_mol_to_ast(input).unwrap();
        CountsValence::new(&counts_model())
            .resolve(&mut ast)
            .unwrap();
        assert_eq!(ast.atoms().count(), atom_count as usize);
        for i in 0..atom_count {
            assert_eq!(ast.atom(AtomId(i)).ast.to_string(), expected_atom);
        }
    }

    #[rstest]
    fn test_parse_smiles_to_ast_methane_counts_resolve() {
        let mut ast = parse_smiles_to_ast("C").unwrap();
        CountsValence::new(&counts_model())
            .resolve(&mut ast)
            .unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.to_string(),
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        );
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
    fn test_parse_mol_bytes_with_resolver_methane_determined() {
        let model = ChemistryModel {
            valence: ValenceModel::Counts(CountsModel {
                table: Cow::Borrowed(ValenceTable::default_table()),
            }),
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
        };
        let ast =
            parse_mol_bytes_with(METHANE_MOL.as_bytes(), &CtfileIoConfig::basic(), &model).unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.to_string(),
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        );
    }

    #[rstest]
    #[case::cfclbri("Br[C@@](F)(Cl)I", 1, vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ])]
    #[case::alanine("N[C@@H](C)C(O)=O", 1, vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
    ])]
    fn test_first_neighbor_toward_ordering(
        #[case] smiles: &str,
        #[case] atom_idx: usize,
        #[case] expected: Vec<StereoLigand>,
    ) {
        let mol = parse_smiles_bytes_to_table_ir(smiles.as_bytes()).unwrap();
        assert_eq!(first_neighbor_toward_ordering(&mol, atom_idx), expected);
    }

    #[rstest]
    #[case::cfclbri("Br[C@@](F)(Cl)I", 1, vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ])]
    #[case::alanine("N[C@@H](C)C(O)=O", 1, vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    ])]
    fn test_tetrahedral_target_ordering(
        #[case] smiles: &str,
        #[case] atom_idx: usize,
        #[case] expected: Vec<StereoLigand>,
    ) {
        let mol = parse_smiles_bytes_to_table_ir(smiles.as_bytes()).unwrap();
        assert_eq!(tetrahedral_target_ordering(&mol, atom_idx), expected);
    }

    #[rstest]
    #[case::cfclbri_clockwise("Br[C@@](F)(Cl)I", 1, 1)]
    #[case::cfclbri_counterclockwise("Br[C@](F)(Cl)I", 1, 0)]
    fn test_raise_tetrahedral_stereo(
        #[case] smiles: &str,
        #[case] atom_idx: usize,
        #[case] expected: u32,
    ) {
        let mol = parse_smiles_bytes_to_table_ir(smiles.as_bytes()).unwrap();
        assert_eq!(
            raise_tetrahedral_stereo(&mol, atom_idx),
            Some(AtomConstraint::TetrahedralStereo(
                StereoConfigurationAst::stereo(StereoCosetAst::Lit(expected))
            ))
        );
    }

    #[rstest]
    #[case::trans("F/C=C/F", 1, 1)]
    #[case::cis("F/C=C\\F", 1, 0)]
    fn test_raise_cis_trans_stereo(
        #[case] smiles: &str,
        #[case] bond_idx: usize,
        #[case] expected: u32,
    ) {
        let mol = parse_smiles_bytes_to_table_ir(smiles.as_bytes()).unwrap();
        assert_eq!(
            raise_cis_trans_stereo(&mol, bond_idx),
            Some(BondConstraint::CisTransStereo(
                StereoConfigurationAst::stereo(StereoCosetAst::Lit(expected))
            ))
        );
    }
}
