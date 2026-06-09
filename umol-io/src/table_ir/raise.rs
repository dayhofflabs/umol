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
    SpinStateAst, StereoConfigurationAst, StereoCosetAst, TryIntoAst, ValueAst,
};
use umol_perm::{space, ClassKey, Permutation};
use umol_shared::error::UmolError;

use crate::table_ir::atom::Atom as TableAtom;
use crate::table_ir::bond::{
    Bond as TableBond, BondDonation as TableBondDonation, BondOrder as TableBondOrder,
};
use crate::table_ir::{BondStereo, Chirality, ChiralityFrame, Molecule as TableMolecule};

mod utils;

use utils::{
    cis_trans_capable, cis_trans_side, first_neighbor_toward_ordering, last_neighbor_away_ordering,
    noncovalent_kind, tetrahedral_ligand_gate, tetrahedral_target_ordering, validate_bond_wedge,
    wedge_winding, wedged_neighbors, StereoBondAtom, StereoFace, StereoLigand,
};

/// Error variants for TableIR -> MoleculeAst raise.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RaiseError {
    /// A tetrahedral stereo assertion on an atom that cannot host a tetrahedral center.
    #[error("tetrahedral stereo at atom {atom} with {count} ligands, expected 3 or 4 ligands")]
    TetrahedralLigandCount { atom: usize, count: usize },
    /// A directional `/`,`\` bond not adjacent to any cis/trans-capable double bond.
    #[error("directional bond {bond} is not adjacent to a stereogenic double bond")]
    DanglingBondWedge { bond: usize },
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

impl TryIntoAst<MoleculeAst> for &TableMolecule {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ast(self, ctx: &Self::Ctx) -> Result<MoleculeAst, RaiseError> {
        let atoms: Vec<AtomAst> = self
            .atoms
            .iter()
            .enumerate()
            .map(|(atom_idx, table_atom)| {
                let mut atom = table_atom.try_into_ast(ctx)?;
                if let Some(constraint) = raise_tetrahedral_stereo(self, atom_idx)? {
                    atom.constraints.add(constraint);
                }
                Ok(atom)
            })
            .collect::<Result<_, _>>()?;

        let mut bonds = Vec::new();
        let mut dative_bonds: Vec<(Vec<AtomId>, AtomId, DativeBondAst)> = Vec::new();
        let mut noncovalent_bonds = Vec::new();
        for (bond_idx, b) in self.bonds.iter().enumerate() {
            validate_bond_wedge(self, bond_idx)?;
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
        let mut atom = AtomAst {
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
        };
        match self.aromatic {
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
        // IO ground defaults for fields the table left unset.
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
        Ok(atom)
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

/// Raise tetrahedral stereo constraint for `atom_idx`.
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

/// Raise cis/trans stereo constraint for `bond_idx`.
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
    // A terminal side (e.g. =O, =CH2) cannot host cis/trans; any flanking marker belongs to a
    // neighboring double bond and is checked by `validate_bond_wedge`.
    if !cis_trans_capable(mol, atom_1, atom_2) {
        return Ok(None);
    }
    let (Some(side_1), Some(side_2)) = (
        cis_trans_side(mol, atom_1, atom_2)?,
        cis_trans_side(mol, atom_2, atom_1)?,
    ) else {
        return Ok(None);
    };
    // Materialize the by-face order from each side's first-ligand face (the second is the opposite).
    let faces = |side: &StereoBondAtom| match side.first_face {
        StereoFace::Above => (side.first_ligand, side.second_ligand),
        StereoFace::Below => (side.second_ligand, side.first_ligand),
    };
    let ((s1_above, s1_below), (s2_above, s2_below)) = (faces(&side_1), faces(&side_2));
    let source = [s1_above, s1_below, s2_above, s2_below];
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
    #[case::ethylideneoxirane_both_ends(parse_smiles_bytes_to_table_ir(b"C/C=C/1CO\\1").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::ethylideneoxirane_open_only(parse_smiles_bytes_to_table_ir(b"C/C=C/1CO1").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::ethylideneoxirane_open_at_oxygen(parse_smiles_bytes_to_table_ir(b"C/C=C(CO\\1)1").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::ethylideneoxirane_open_at_oxygen_both_ends(parse_smiles_bytes_to_table_ir(b"C/C=C(CO\\1)/1").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::ethylideneoxirane_open_at_oxygen_close_only(parse_smiles_bytes_to_table_ir(b"C/C=C(CO1)/1").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::cyclooctene_trans(parse_smiles_bytes_to_table_ir(b"C1=C/CCCCCC/1").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::cyclooctene_trans_open_marked(parse_smiles_bytes_to_table_ir(b"C\\1=C/CCCCCC1").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::cyclooctene_cis(parse_smiles_bytes_to_table_ir(b"C1=C\\CCCCCC/1").unwrap(), 1, Some(StereoCosetAst::Lit(0)))]
    #[case::geminal_difluoro(parse_smiles_bytes_to_table_ir(b"F/C(F)=C(C)\\CC").unwrap(), 2, Some(StereoCosetAst::Lit(1)))]
    #[case::butanone_oxime(parse_smiles_bytes_to_table_ir(b"C/C(CC)=N\\O").unwrap(), 3, Some(StereoCosetAst::Lit(0)))]
    #[case::fluoropropene_e_backslash(parse_smiles_bytes_to_table_ir(b"F\\C=C\\C").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::fluoropropene_e_methyl_first(parse_smiles_bytes_to_table_ir(b"C/C=C/F").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::fluoropropene_e_methyl_first_backslash(parse_smiles_bytes_to_table_ir(b"C\\C=C\\F").unwrap(), 1, Some(StereoCosetAst::Lit(1)))]
    #[case::trisubstituted(parse_smiles_bytes_to_table_ir(b"F/C(C)=C(Cl)/C").unwrap(), 2, Some(StereoCosetAst::Lit(0)))]
    #[case::mol_either(parse_mol_bytes_to_table_ir(CIS_TRANS_EITHER_MOL.as_bytes()).unwrap(), 1, Some(StereoCosetAst::Undetermined))]
    #[case::one_sided_marker(parse_smiles_bytes_to_table_ir(b"C(C)=C(Cl)/C").unwrap(), 1, None)]
    #[case::plain_double(parse_smiles_bytes_to_table_ir(b"C=C").unwrap(), 0, None)]
    #[case::terminal_no_substituent(parse_smiles_bytes_to_table_ir(b"F/C=C").unwrap(), 1, None)]
    #[case::cyclohexenone_carbonyl(parse_smiles_bytes_to_table_ir(b"O=C1/C=C\\CCC1").unwrap(), 0, None)]
    #[case::cyclohexenone(parse_smiles_bytes_to_table_ir(b"O=C1/C=C\\CCC1").unwrap(), 3, Some(StereoCosetAst::Lit(0)))]
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
    fn test_raise_cis_trans_stereo_error(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: RaiseError,
    ) {
        assert_eq!(raise_cis_trans_stereo(&mol, bond_idx), Err(expected));
    }

    #[rstest]
    #[case::dangling(parse_smiles_bytes_to_table_ir(b"F/C=C").unwrap(), 0, Err(RaiseError::DanglingBondWedge { bond: 0 }))]
    #[case::flanks_capable(parse_smiles_bytes_to_table_ir(b"O=C1/C=C\\CCC1").unwrap(), 2, Ok(()))]
    fn test_validate_bond_wedge(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: Result<(), RaiseError>,
    ) {
        assert_eq!(validate_bond_wedge(&mol, bond_idx), expected);
    }
}
