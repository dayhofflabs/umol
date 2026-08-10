//! TableIR → `umol_graph_ir::Molecule` raise.
//!
//! Implements `TryIntoIr<Molecule> for &Molecule` (and the per-atom and
//! per-bond analogues). Table IR fields copy to `Lit` / `Undetermined`; IO
//! raise applies fixed IO ground semantics for resolution.

use std::any::Any;
use std::collections::HashSet;

use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId, BondConstraintForm, BondForm,
    BooleanForm, CisTransStereoForm, Constraints, DativeBondForm, ElementForm, IsotopeMassForm,
    Lattice, Molecule, MoleculeEntries, MoleculeEntriesError, MulticenterBondForm,
    NoncovalentBondForm, NumForm, StereoCoset, TetrahedralStereoForm, TryIntoIr,
    UnpairedElectronsForm,
};
use umol_perm::{ClassKey, Permutation};
use umol_utils::error::UmolError;

use crate::table_ir::atom::Atom as TableAtom;
use crate::table_ir::bond::{
    Bond as TableBond, BondDonation as TableBondDonation, BondOrder as TableBondOrder,
};
use crate::table_ir::raise::utils::coset_from_wedge_winding;
use crate::table_ir::{BondStereo, Chirality, ChiralityFrame, Molecule as TableMolecule};

mod utils;

use utils::{
    cis_trans_capable, cis_trans_side, first_neighbor_toward_ordering, last_neighbor_away_ordering,
    neighbor_count, noncovalent_kind, tetrahedral_ligand_ordering, validate_bond_direction,
    validate_tetrahedral_geometry, wedge_bond_neighbors, StereoBondAtom, StereoHalfplane,
};

/// Error variants for TableIR -> Molecule raise.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RaiseError {
    #[error(transparent)]
    MoleculeEntries(#[from] MoleculeEntriesError),
    #[error("tetrahedral stereo at atom {atom} with {count} ligands, expected 3 or 4 ligands")]
    TetrahedralLigandCount { atom: usize, count: usize },
    #[error("directional bond {bond} not adjacent to a stereogenic double bond")]
    DanglingBondDirection { bond: usize },
    #[error("contradictory cis/trans markers at atom {atom}")]
    CisTransConflict { atom: usize },
    #[error("inconsistent wedge bonds at atom {atom}")]
    WedgeConflict { atom: usize },
}

impl UmolError for RaiseError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl TryIntoIr<Molecule> for &TableMolecule {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ir(self, ctx: &Self::Ctx) -> Result<Molecule, RaiseError> {
        let atoms: Vec<AtomForm> = self
            .atoms
            .iter()
            .enumerate()
            .map(|(atom_idx, table_atom)| {
                let mut atom = table_atom.try_into_ir(ctx)?;
                if let Some(constraint) = raise_tetrahedral_stereo(self, atom_idx)? {
                    atom.constraints.set(constraint);
                }
                Ok(atom)
            })
            .collect::<Result<_, RaiseError>>()?;

        let mut bonds = Vec::new();
        let mut dative_bonds: Vec<(Vec<AtomId>, AtomId, DativeBondForm)> = Vec::new();
        let mut noncovalent_bonds = Vec::new();
        for (bond_idx, b) in self.bonds.iter().enumerate() {
            validate_bond_direction(self, bond_idx)?;
            let a_idx = AtomId(b.atoms.first());
            let b_idx = AtomId(b.atoms.second());
            if let Some(kind) = b.noncovalent.map(noncovalent_kind) {
                noncovalent_bonds.push((a_idx, b_idx, NoncovalentBondForm::from_kind(kind)));
            } else if let Some(donation) = b.donation {
                let (donor, acceptor) = match donation {
                    TableBondDonation::Donating => (a_idx, b_idx),
                    TableBondDonation::Accepting => (b_idx, a_idx),
                    _ => {
                        bonds.push((a_idx, b_idx, b.try_into_ir(ctx)?));
                        continue;
                    }
                };
                let dative_bond = DativeBondForm::new(raise_bond_order(b.order));
                dative_bonds.push((vec![donor], acceptor, dative_bond));
            } else {
                let mut bond_form = b.try_into_ir(ctx)?;
                if let Some(constraint) = raise_cis_trans_stereo(self, bond_idx)? {
                    bond_form.constraints.set(constraint);
                }
                bonds.push((a_idx, b_idx, bond_form));
            }
        }

        let multicenter_bond: Vec<(Vec<AtomId>, MulticenterBondForm)> = self
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
                (atoms, MulticenterBondForm::default())
            })
            .collect();

        let constraints = Constraints::new();

        Molecule::try_from_entries(MoleculeEntries {
            atoms,
            bonds,
            dative: dative_bonds,
            multicenter: multicenter_bond,
            noncovalent: noncovalent_bonds,
            constraints,
            ..Default::default()
        })
        .map_err(Into::into)
    }
}

impl TryIntoIr<AtomForm> for &TableAtom {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ir(self, _ctx: &Self::Ctx) -> Result<AtomForm, RaiseError> {
        let mut atom = AtomForm {
            element: match self.element {
                Some(element) => ElementForm::Lit(element),
                None => ElementForm::Undetermined,
            },
            isotope_mass: match self.isotope_mass {
                Some(m) => IsotopeMassForm::Lit(m),
                None => IsotopeMassForm::Undetermined,
            },
            charge: match self.charge {
                Some(c) => NumForm::Lit(c as i64),
                None => NumForm::Undetermined,
            },
            implicit_hydrogens: match self.implicit_hydrogens {
                Some(n) => NumForm::Lit(n as i64),
                None => NumForm::Undetermined,
            },
            lone_pairs: match self.lone_pairs {
                Some(n) => NumForm::Lit(n as i64),
                None => NumForm::Undetermined,
            },
            unpaired_electrons: UnpairedElectronsForm {
                count: match self.unpaired_electrons {
                    Some(unpaired_electrons) => NumForm::Lit(unpaired_electrons as i64),
                    None => NumForm::Undetermined,
                },
                multiplicity: match self.multiplicity {
                    Some(multiplicity) => NumForm::Lit(u8::from(multiplicity) as i64),
                    None => NumForm::Undetermined,
                },
            },
            constraints: Default::default(),
        };
        match self.aromatic {
            Some(true) => {
                atom.constraints.set(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(NumForm::Undetermined),
                ));
                // A bare aromatic heteroatom specifies zero H; any H must be bracketed
                // ([nH]), which arrives above as an explicit count.
                if matches!(self.element, Some(element) if element != Element::C)
                    && matches!(atom.implicit_hydrogens, NumForm::Undetermined)
                {
                    atom.implicit_hydrogens = NumForm::Lit(0);
                }
            }
            Some(false) => {
                atom.constraints.set(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::NotAromatic,
                ));
            }
            None => {}
        }
        // IO ground defaults for fields the table left unset.
        if matches!(atom.isotope_mass, IsotopeMassForm::Undetermined) {
            atom.isotope_mass = IsotopeMassForm::Natural;
        }
        if matches!(atom.charge, NumForm::Undetermined) {
            atom.charge = NumForm::Lit(0);
        }
        if matches!(atom.unpaired_electrons.count, NumForm::Undetermined) {
            atom.unpaired_electrons.count = NumForm::Lit(0);
        }
        atom.constraints.retain(|c| !c.is_undetermined());
        Ok(atom)
    }
}

impl TryIntoIr<BondForm> for &TableBond {
    type Ctx = ();
    type Error = RaiseError;

    fn try_into_ir(self, _ctx: &Self::Ctx) -> Result<BondForm, RaiseError> {
        let mut bond = BondForm::new(raise_bond_order(self.order));
        bond.charge = match self.charge {
            Some(c) => NumForm::Lit(c as i64),
            None => NumForm::Undetermined,
        };
        if matches!(self.order, TableBondOrder::Aromatic) {
            bond.constraints
                .set(BondConstraintForm::Aromatic(BooleanForm::Lit(true)));
        }
        Ok(bond)
    }
}

fn raise_bond_order(order: TableBondOrder) -> NumForm {
    match order {
        TableBondOrder::Zero => NumForm::Lit(0),
        TableBondOrder::Single => NumForm::Lit(1),
        TableBondOrder::Double => NumForm::Lit(2),
        TableBondOrder::Triple => NumForm::Lit(3),
        TableBondOrder::Quadruple => NumForm::Lit(4),
        TableBondOrder::Quintuple => NumForm::Lit(5),
        TableBondOrder::Sextuple => NumForm::Lit(6),
        // Definite-aromatic: localized bond order is 1 by Kekulé convention;
        // the aromatic flag is added separately as `BondConstraintForm::Aromatic`.
        // Renders as `1#a`.
        TableBondOrder::Aromatic => NumForm::Lit(1),
        // Fuzzy orders: no concrete bond order can be assigned; raise to
        // `Undetermined`. Aromatic-flag setting (where applicable) is left
        // off — the chemistry of these is too ambiguous for the raise.
        // TODO: Convert to LitSet.
        TableBondOrder::SingleOrDouble
        | TableBondOrder::SingleOrAromatic
        | TableBondOrder::DoubleOrAromatic
        | TableBondOrder::Any => NumForm::Undetermined,
    }
}

/// Raise tetrahedral stereo constraint for `atom_idx`.
fn raise_tetrahedral_stereo(
    mol: &TableMolecule,
    atom_idx: usize,
) -> Result<Option<AtomConstraintForm>, RaiseError> {
    let chirality = mol.atoms[atom_idx].chirality;
    let (relabeling, source_coset): (Permutation, usize) = match chirality {
        Some(Chirality::Unspecified) => {
            validate_tetrahedral_geometry(mol, atom_idx)?;
            return Ok(Some(AtomConstraintForm::TetrahedralStereo(
                TetrahedralStereoForm::stereo(StereoCoset::Undetermined),
            )));
        }
        Some(symbol) => {
            let source_coset = match (symbol, mol.chirality_frame) {
                (
                    Chirality::CounterClockwise | Chirality::Tetrahedral { arr: 1 },
                    Some(ChiralityFrame::FirstNeighborToward),
                ) => 0,
                (
                    Chirality::Clockwise | Chirality::Tetrahedral { arr: 2 },
                    Some(ChiralityFrame::FirstNeighborToward),
                ) => 1,
                (
                    Chirality::CounterClockwise | Chirality::Tetrahedral { arr: 1 },
                    Some(ChiralityFrame::LastNeighborAway),
                ) => 1,
                (
                    Chirality::Clockwise | Chirality::Tetrahedral { arr: 2 },
                    Some(ChiralityFrame::LastNeighborAway),
                ) => 0,
                _ => return Ok(None),
            };
            validate_tetrahedral_geometry(mol, atom_idx)?;
            let source_ordering = match mol.chirality_frame {
                Some(ChiralityFrame::FirstNeighborToward) => {
                    first_neighbor_toward_ordering(mol, atom_idx)
                }
                Some(ChiralityFrame::LastNeighborAway) => {
                    last_neighbor_away_ordering(mol, atom_idx)
                }
                None => return Ok(None),
            };
            let permutation = Permutation::between(
                &source_ordering,
                &tetrahedral_ligand_ordering(mol, atom_idx),
            )
            .expect("validated tetrahedral frames contain the same ligands");
            (permutation, source_coset)
        }
        None => {
            let Some(positions) = mol.positions.as_ref() else {
                return Ok(None);
            };
            // Exclude atoms adjacent to tetrahedral stereo centers that share a wedge bond.
            let count = neighbor_count(mol, atom_idx);
            if count != 3 && count != 4 {
                return Ok(None);
            }
            let neighbors = wedge_bond_neighbors(mol, atom_idx);
            let Some(&(neighbor_idx, outofplane)) = neighbors.first() else {
                return Ok(None);
            };
            let target_ordering = tetrahedral_ligand_ordering(mol, atom_idx);
            let source_coset = coset_from_wedge_winding(
                atom_idx,
                &target_ordering,
                neighbor_idx,
                positions,
                outofplane,
            );
            for &(neighbor_idx, outofplane) in &neighbors[1..] {
                if coset_from_wedge_winding(
                    atom_idx,
                    &target_ordering,
                    neighbor_idx,
                    positions,
                    outofplane,
                ) != source_coset
                {
                    return Err(RaiseError::WedgeConflict { atom: atom_idx });
                }
            }
            (Permutation::identity(4), source_coset)
        }
    };
    let coset = ClassKey::Tetrahedral
        .space()
        .reindex(source_coset as u32, relabeling)
        .expect("tetrahedral coset reindex");
    Ok(Some(AtomConstraintForm::TetrahedralStereo(
        TetrahedralStereoForm::stereo(StereoCoset::Lit(coset)),
    )))
}

/// Raise cis/trans stereo constraint for `bond_idx`.
fn raise_cis_trans_stereo(
    mol: &TableMolecule,
    bond_idx: usize,
) -> Result<Option<BondConstraintForm>, RaiseError> {
    let bond = &mol.bonds[bond_idx];
    if bond.order != TableBondOrder::Double {
        return Ok(None);
    }
    if bond.stereo == Some(BondStereo::Either) {
        return Ok(Some(BondConstraintForm::CisTransStereo(
            CisTransStereoForm::stereo(StereoCoset::Undetermined),
        )));
    }
    let atom_1_idx = bond.start_atom() as usize;
    let atom_2_idx = bond.end_atom() as usize;
    if !cis_trans_capable(mol, atom_1_idx, atom_2_idx) {
        return Ok(None);
    }
    let (Some(side_1), Some(side_2)) = (
        cis_trans_side(mol, atom_1_idx, atom_2_idx)?,
        cis_trans_side(mol, atom_2_idx, atom_1_idx)?,
    ) else {
        return Ok(None);
    };
    // Generate the halfplane assignments for each side of the double bond.
    let halfplanes = |side: &StereoBondAtom| match side.first_halfplane {
        StereoHalfplane::Top => (side.first_ligand, side.second_ligand),
        StereoHalfplane::Bottom => (side.second_ligand, side.first_ligand),
    };
    let ((s1_above, s1_below), (s2_above, s2_below)) = (halfplanes(&side_1), halfplanes(&side_2));
    let source = [s1_above, s1_below, s2_above, s2_below];
    let target = [
        side_1.first_ligand,
        side_1.second_ligand,
        side_2.first_ligand,
        side_2.second_ligand,
    ];
    let coset = ClassKey::CisTrans
        .space()
        .index(
            Permutation::between(&source, &target)
                .expect("validated cis/trans frames contain the same ligands"),
        )
        .expect("cis/trans coset index");
    Ok(Some(BondConstraintForm::CisTransStereo(
        CisTransStereoForm::stereo(StereoCoset::Lit(coset)),
    )))
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_chem::spin::SpinMultiplicity;
    use umol_graph_ir::ir::{AtomConstraintsForm, BondId, Entity};

    use super::*;
    use crate::ctfile::parse_mol_to_ir;
    use crate::ctfile::parser::parse_mol_bytes_to_table_ir;
    use crate::smiles::Smiles;
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

    const WEDGE_CONFLICT_MOL: &str = "wedge-conflict\n\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  1  0  0  0\n  1  3  1  1  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n";

    const CFCLBRI_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.6906   -0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.6906    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.6906    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  6        0\n  2  4  1  0        0\n  2  5  1  6        0\nM  END\n";

    const CFCLBRI_R_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.6906   -0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.6906    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.6906    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  1        0\n  2  4  1  0        0\n  2  5  1  1        0\nM  END\n";

    const CFCLBRI_INCONSISTENT_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    0.6906   -0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.6906    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.6906    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  6        0\n  2  4  1  0        0\n  2  5  1  1        0\nM  END\n";

    const CFCLBRI_SINGLE_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.6906   -0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.6906    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.6906    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  0        0\n  2  4  1  0        0\n  2  5  1  1        0\nM  END\n";

    // Gives opposite coset from SMILES `C[C@@H](O)CC` example because of atom ordering swap.
    const BUTANOL_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n   -1.0643   -0.6145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3548   -0.2048    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3548   -0.6145    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3548    0.6145    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0643   -0.2048    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  0        0\n  2  4  1  1        0\n  3  5  1  0        0\nM  END\n";

    // Same coset as `C[C@H](N)C(O)=O` because atom ordering differs by even (3-cycle) permutation.
    const ALANINE_WEDGE_MOL: &str = "\n\n\n  6  5  0  0  1  0  0  0  0  0999 V2000\n   -0.3560   -1.0277    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0680    0.2055    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3560    1.0277    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3560   -0.2055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.3560    0.2055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0680   -0.2055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  4  1  2  0        0\n  4  2  1  0        0\n  5  3  1  1        0\n  5  4  1  0        0\n  6  5  1  0        0\nM  END\n";

    // Opposite coset from `C[S@@+]([O-])CC` because of atom ordering swap.
    const SULFOXIDE_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n   -1.0680   -0.6166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3560   -0.2055    0.0000 S   0  3  0  0  0  0  0  0  0  0  0  0\n    0.3560   -0.6166    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3560    0.6166    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n    1.0680   -0.2055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  0        0\n  2  4  1  1        0\n  3  5  1  0        0\nM  CHG  2   2   1   4  -1\nM  END\n";

    // 2-methyloxirane, explicit H, two consistent wedges. Opposite coset from `C[C@@]1([H])OC1`.
    const METHYLOXIRANE_WEDGE_MOL: &str = "\n\n\n  5  5  0  0  1  0  0  0  0  0999 V2000\n   -0.1738    0.0355    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4076    0.6168    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9889    0.1428    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.9889    0.0355    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6743   -0.6168    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  1  3  1  1        0\n  2  4  1  0        0\n  4  1  1  0        0\n  1  5  1  6        0\nM  END\n";

    const PROCHIRAL_METHYLENE_WEDGE_MOL: &str = "\n\n\n  7  6  0  0  0  0  0  0  0  0999 V2000\n   -0.3009   -0.2055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4111    0.2055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4111    1.0277    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.1231   -0.2055    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.3009   -1.0277    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7120    0.5065    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.1231   -0.2055    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  0        0\n  2  4  2  0        0\n  1  5  1  0        0\n  1  6  1  1        0\n  1  7  1  6        0\nM  END\n";

    const CHIRAL_PARITY_MOL: &str = "chiral\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n";

    const CIS_TRANS_EITHER_MOL: &str = "butene\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  3  0  0  0\n  3  4  1  0  0  0  0\nM  END\n";

    #[rstest]
    fn test_table_molecule_try_into_ir(methane: TableMolecule) {
        let molecule: Molecule = (&methane).try_into_ir(&()).unwrap();
        assert_eq!(
            molecule,
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm {
                    element: ElementForm::Lit(Element::C),
                    isotope_mass: IsotopeMassForm::Natural,
                    charge: NumForm::Lit(0),
                    implicit_hydrogens: NumForm::Lit(4),
                    lone_pairs: NumForm::Undetermined,
                    unpaired_electrons: UnpairedElectronsForm {
                        count: NumForm::Lit(0),
                        multiplicity: NumForm::Undetermined,
                    },
                    constraints: AtomConstraintsForm::new(),
                }],
                ..Default::default()
            })
        );
    }

    #[rstest]
    #[case::shared_cis_trans_ligand(
        Smiles::parse_bytes(b"SSC=S1CC1\\2C=112").unwrap().into_table_ir(),
        RaiseError::DanglingBondDirection { bond: 6 }
    )]
    #[case::invalid_bond_endpoint(
        {
            let mut molecule = TableMolecule::empty();
            molecule.atoms.push(TableAtom::from_element(Element::C));
            molecule.bonds.push(TableBond::new(0, 1, TableBondOrder::Single));
            molecule
        },
        RaiseError::MoleculeEntries(MoleculeEntriesError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        })
    )]
    fn test_table_molecule_try_into_ir_error(
        #[case] molecule: TableMolecule,
        #[case] expected: RaiseError,
    ) {
        let actual: Result<Molecule, RaiseError> = (&molecule).try_into_ir(&());
        assert_eq!(actual, Err(expected));
    }

    #[rstest]
    #[case::bare(
        None,
        None,
        None,
        None,
        None,
        None,
        AtomForm {
            element: ElementForm::Undetermined,
            isotope_mass: IsotopeMassForm::Natural,
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Undetermined,
            lone_pairs: NumForm::Undetermined,
            unpaired_electrons: UnpairedElectronsForm {
                count: NumForm::Lit(0),
                multiplicity: NumForm::Undetermined,
            },
            constraints: AtomConstraintsForm::new(),
        }
    )]
    #[case::bracket_fields(
        Some(13),
        Some(-1),
        Some(2),
        Some(1),
        Some(2),
        Some(SpinMultiplicity::SINGLET),
        AtomForm {
            element: ElementForm::Undetermined,
            isotope_mass: IsotopeMassForm::Lit(13),
            charge: NumForm::Lit(-1),
            implicit_hydrogens: NumForm::Lit(2),
            lone_pairs: NumForm::Lit(1),
            unpaired_electrons: UnpairedElectronsForm {
                count: NumForm::Lit(2),
                multiplicity: NumForm::Lit(1),
            },
            constraints: AtomConstraintsForm::new(),
        }
    )]
    fn test_table_atom_try_into_ir(
        #[case] isotope_mass: Option<u32>,
        #[case] charge: Option<i8>,
        #[case] implicit_hydrogens: Option<u8>,
        #[case] lone_pairs: Option<u8>,
        #[case] unpaired_electrons: Option<u8>,
        #[case] multiplicity: Option<SpinMultiplicity>,
        #[case] expected: AtomForm,
    ) {
        let atom = TableAtom {
            isotope_mass,
            charge,
            implicit_hydrogens,
            lone_pairs,
            unpaired_electrons,
            multiplicity,
            ..TableAtom::wildcard()
        };
        assert_eq!(atom.try_into_ir(&()), Ok(expected));
    }

    #[rstest]
    fn test_table_atom_try_into_ir_aromatic_wildcard() {
        let atom = TableAtom {
            aromatic: Some(true),
            ..TableAtom::wildcard()
        };
        assert_eq!(
            atom.try_into_ir(&()),
            Ok(AtomForm {
                element: ElementForm::Undetermined,
                isotope_mass: IsotopeMassForm::Natural,
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Undetermined,
                lone_pairs: NumForm::Undetermined,
                unpaired_electrons: UnpairedElectronsForm {
                    count: NumForm::Lit(0),
                    multiplicity: NumForm::Undetermined,
                },
                constraints: AtomConstraintsForm::from(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(NumForm::Undetermined),
                )),
            })
        );
    }

    #[rstest]
    #[case::table_aromatic_none(None, None)]
    #[case::table_aromatic_false(Some(false), Some(AromaticValenceForm::NotAromatic))]
    #[case::table_aromatic_true(
        Some(true),
        Some(AromaticValenceForm::Aromatic(NumForm::Undetermined))
    )]
    fn test_table_molecule_try_into_ir_aromatic(
        mut carbon: TableMolecule,
        #[case] aromatic: Option<bool>,
        #[case] expected: Option<AromaticValenceForm>,
    ) {
        carbon.atoms[0].aromatic = aromatic;
        let molecule: Molecule = (&carbon).try_into_ir(&()).unwrap();
        assert_eq!(
            molecule
                .atom(AtomId(0))
                .attributes
                .constraints
                .aromatic_valence(),
            expected.as_ref()
        );
    }

    // A bare aromatic heteroatom resolves to zero H; aromatic carbon and explicit
    // bracket H are left to the valence model / preserved.
    #[rstest]
    #[case::aromatic_nitrogen_bare(Element::N, Some(true), None, NumForm::Lit(0))]
    #[case::aromatic_oxygen_bare(Element::O, Some(true), None, NumForm::Lit(0))]
    #[case::aromatic_nitrogen_bracket_h(Element::N, Some(true), Some(1), NumForm::Lit(1))]
    #[case::aromatic_carbon_bare(Element::C, Some(true), None, NumForm::Undetermined)]
    #[case::aliphatic_nitrogen_bare(Element::N, Some(false), None, NumForm::Undetermined)]
    fn test_table_molecule_try_into_ir_aromatic_heteroatoms(
        #[case] element: Element,
        #[case] aromatic: Option<bool>,
        #[case] hydrogens: Option<u8>,
        #[case] expected: NumForm,
    ) {
        let mut atom = TableAtom::from_element(element);
        atom.aromatic = aromatic;
        atom.implicit_hydrogens = hydrogens;
        let mut mol = TableMolecule::empty();
        mol.atoms.push(atom);
        let molecule: Molecule = (&mol).try_into_ir(&()).unwrap();
        assert_eq!(
            molecule.atom(AtomId(0)).attributes.implicit_hydrogens,
            expected
        );
    }

    #[rstest]
    fn test_table_molecule_try_into_ir_bond_order(
        #[with(TableBondOrder::Double)] diatomic: TableMolecule,
    ) {
        let molecule: Molecule = (&diatomic).try_into_ir(&()).unwrap();
        let bond = molecule.bond(BondId(0)).attributes;
        assert!(matches!(bond.order, NumForm::Lit(2)));
    }

    #[rstest]
    fn test_table_molecule_try_into_ir_aromatic_bond(
        #[with(TableBondOrder::Aromatic)] diatomic: TableMolecule,
    ) {
        let molecule: Molecule = (&diatomic).try_into_ir(&()).unwrap();
        let bond = molecule.bond(BondId(0)).attributes;
        assert!(matches!(bond.order, NumForm::Lit(1)));
        assert!(bond
            .constraints
            .iter()
            .any(|c| matches!(c, BondConstraintForm::Aromatic(BooleanForm::Lit(true)))));
        for i in 0..2 {
            assert!(molecule
                .atom(AtomId(i))
                .attributes
                .constraints
                .aromatic_valence()
                .is_none());
        }
    }

    #[rstest]
    #[case::methane(METHANE_MOL, "C#i=#c0#u0")]
    #[case::benzene(BENZENE_AROMATIC_MOL, "C#i=#c0#u0")]
    #[case::carbon_h0(CARBON_H0_EXPLICIT_MOL, "C#i=#c0#h0#u0")]
    fn test_parse_mol_to_ir(#[case] input: &str, #[case] expected_atom: &str) {
        let molecule = parse_mol_to_ir(input).unwrap();
        let atom = molecule.atom(AtomId(0)).attributes;
        assert_eq!(atom.charge, NumForm::Lit(0));
        assert!(atom.constraints.aromatic_valence().is_none());
        assert_eq!(atom.to_string(), expected_atom);
    }

    #[rstest]
    #[case::organic("C", "C#i=#c0#u0#a!")]
    fn test_table_molecule_try_into_ir_smiles(#[case] input: &str, #[case] expected_atom: &str) {
        let smiles = Smiles::parse(input).unwrap();
        let molecule: Molecule = smiles.as_table_ir().try_into_ir(&()).unwrap();
        let atom = molecule.atom(AtomId(0)).attributes;
        assert_eq!(atom.charge, NumForm::Lit(0));
        assert!(matches!(atom.implicit_hydrogens, NumForm::Undetermined));
        assert!(matches!(
            atom.constraints.aromatic_valence(),
            Some(AromaticValenceForm::NotAromatic)
        ));
        assert_eq!(atom.to_string(), expected_atom);
    }

    #[rstest]
    fn test_table_molecule_try_into_ir_smiles_wildcard() {
        let smiles = Smiles::parse("*").unwrap();
        let molecule: Molecule = smiles.as_table_ir().try_into_ir(&()).unwrap();

        assert_eq!(
            molecule.atom(AtomId(0)).attributes,
            &AtomForm {
                element: ElementForm::Undetermined,
                isotope_mass: IsotopeMassForm::Natural,
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Undetermined,
                lone_pairs: NumForm::Undetermined,
                unpaired_electrons: UnpairedElectronsForm {
                    count: NumForm::Lit(0),
                    multiplicity: NumForm::Undetermined,
                },
                constraints: AtomConstraintsForm::new(),
            }
        );
    }

    #[rstest]
    #[case::cfclbri_clockwise(Smiles::parse_bytes(b"Br[C@@](F)(Cl)I").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::cfclbri_counterclockwise(Smiles::parse_bytes(b"Br[C@](F)(Cl)I").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::cfclbri_fluorine_first(Smiles::parse_bytes(b"F[C@](Cl)(Br)I").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::methyloxirane_explicit_h(Smiles::parse_bytes(b"C[C@@]1([H])OC1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::butan_2_ol(Smiles::parse_bytes(b"C[C@@H](O)CC").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::alanine(Smiles::parse_bytes(b"C[C@H](N)C(O)=O").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::ring_then_branch(Smiles::parse_bytes(b"C[C@]1(Cl)CC(C)CC1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::branch_then_ring(Smiles::parse_bytes(b"C[C@](Cl)1CC(C)CC1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_parity_clockwise(parse_mol_bytes_to_table_ir(CHIRAL_PARITY_MOL.as_bytes()).unwrap(), 0, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_cfclbri(parse_mol_bytes_to_table_ir(CFCLBRI_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_cfclbri_r(parse_mol_bytes_to_table_ir(CFCLBRI_R_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_wedge_cfclbri_single(parse_mol_bytes_to_table_ir(CFCLBRI_SINGLE_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_wedge_butanol(parse_mol_bytes_to_table_ir(BUTANOL_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_alanine(parse_mol_bytes_to_table_ir(ALANINE_WEDGE_MOL.as_bytes()).unwrap(), 4, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_sulfoxide(parse_mol_bytes_to_table_ir(SULFOXIDE_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_methyloxirane(parse_mol_bytes_to_table_ir(METHYLOXIRANE_WEDGE_MOL.as_bytes()).unwrap(), 0, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_prochiral_methylene(parse_mol_bytes_to_table_ir(PROCHIRAL_METHYLENE_WEDGE_MOL.as_bytes()).unwrap(), 0, Some(StereoCoset::Lit(1)))]
    #[case::sulfoxide_counterclockwise(Smiles::parse_bytes(b"C[S@](=O)CC").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::sulfoxide_clockwise(Smiles::parse_bytes(b"C[S@@](=O)CC").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::sulfoxide_charge_separated(Smiles::parse_bytes(b"C[S@@+]([O-])CC").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::no_descriptor(Smiles::parse_bytes(b"F[C@](Cl)(Br)I").unwrap().into_table_ir(), 0, None)]
    fn test_raise_tetrahedral_stereo(
        #[case] mol: TableMolecule,
        #[case] atom_idx: usize,
        #[case] expected: Option<StereoCoset>,
    ) {
        let expected = expected.map(|coset| {
            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::stereo(coset))
        });
        assert_eq!(raise_tetrahedral_stereo(&mol, atom_idx), Ok(expected));
    }

    #[rstest]
    #[case::dimethyl_sulfide(Smiles::parse_bytes(b"C[S@]C").unwrap().into_table_ir(), 1, RaiseError::TetrahedralLigandCount { atom: 1, count: 2 })]
    #[case::parallel_bonds(Smiles::parse_bytes(b"C[C]2[C@@]2[C-]").unwrap().into_table_ir(), 2, RaiseError::TetrahedralLigandCount { atom: 2, count: 2 })]
    #[case::wedge_conflict(parse_mol_bytes_to_table_ir(WEDGE_CONFLICT_MOL.as_bytes()).unwrap(), 0, RaiseError::WedgeConflict { atom: 0 })]
    #[case::cfclbri_inconsistent_wedges(parse_mol_bytes_to_table_ir(CFCLBRI_INCONSISTENT_WEDGE_MOL.as_bytes()).unwrap(), 1, RaiseError::WedgeConflict { atom: 1 })]
    fn test_raise_tetrahedral_stereo_error(
        #[case] mol: TableMolecule,
        #[case] atom_idx: usize,
        #[case] expected: RaiseError,
    ) {
        assert_eq!(raise_tetrahedral_stereo(&mol, atom_idx), Err(expected));
    }

    #[rstest]
    #[case::trans(Smiles::parse_bytes(b"F/C=C/F").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::cis(Smiles::parse_bytes(b"F/C=C\\F").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::fluoropropene_e(Smiles::parse_bytes(b"F/C=C/C").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::fluoropropene_z(Smiles::parse_bytes(b"F/C=C\\C").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::fluoropropene_z_flipped(Smiles::parse_bytes(b"F\\C=C/C").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::fluoropropene_z_methyl_first(Smiles::parse_bytes(b"C/C=C\\F").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::azomethane_e(Smiles::parse_bytes(b"C/N=N/C").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::azomethane_z(Smiles::parse_bytes(b"C/N=N\\C").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::ethylideneoxirane(Smiles::parse_bytes(b"C/C=C1CO\\1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::ethylideneoxirane_both_ends(Smiles::parse_bytes(b"C/C=C/1CO\\1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::ethylideneoxirane_open_only(Smiles::parse_bytes(b"C/C=C/1CO1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::ethylideneoxirane_open_at_oxygen(Smiles::parse_bytes(b"C/C=C(CO\\1)1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::ethylideneoxirane_open_at_oxygen_both_ends(Smiles::parse_bytes(b"C/C=C(CO\\1)/1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::ethylideneoxirane_open_at_oxygen_close_only(Smiles::parse_bytes(b"C/C=C(CO1)/1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::cyclooctene_trans(Smiles::parse_bytes(b"C1=C/CCCCCC/1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::cyclooctene_trans_open_marked(Smiles::parse_bytes(b"C\\1=C/CCCCCC1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::cyclooctene_cis(Smiles::parse_bytes(b"C1=C\\CCCCCC/1").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(0)))]
    #[case::geminal_difluoro(Smiles::parse_bytes(b"F/C(F)=C(C)\\CC").unwrap().into_table_ir(), 2, Some(StereoCoset::Lit(1)))]
    #[case::butanone_oxime(Smiles::parse_bytes(b"C/C(CC)=N\\O").unwrap().into_table_ir(), 3, Some(StereoCoset::Lit(0)))]
    #[case::fluoropropene_e_backslash(Smiles::parse_bytes(b"F\\C=C\\C").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::fluoropropene_e_methyl_first(Smiles::parse_bytes(b"C/C=C/F").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::fluoropropene_e_methyl_first_backslash(Smiles::parse_bytes(b"C\\C=C\\F").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::trisubstituted(Smiles::parse_bytes(b"F/C(C)=C(Cl)/C").unwrap().into_table_ir(), 2, Some(StereoCoset::Lit(0)))]
    #[case::mol_either(parse_mol_bytes_to_table_ir(CIS_TRANS_EITHER_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Undetermined))]
    #[case::one_sided_marker(Smiles::parse_bytes(b"C(C)=C(Cl)/C").unwrap().into_table_ir(), 1, None)]
    #[case::plain_double(Smiles::parse_bytes(b"C=C").unwrap().into_table_ir(), 0, None)]
    #[case::terminal_no_substituent(Smiles::parse_bytes(b"F/C=C").unwrap().into_table_ir(), 1, None)]
    #[case::cyclohexenone_carbonyl(Smiles::parse_bytes(b"O=C1/C=C\\CCC1").unwrap().into_table_ir(), 0, None)]
    #[case::cyclohexenone(Smiles::parse_bytes(b"O=C1/C=C\\CCC1").unwrap().into_table_ir(), 3, Some(StereoCoset::Lit(0)))]
    #[case::hexadiene_first(Smiles::parse_bytes(b"C/C=C/C=C/C").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::hexadiene_second(Smiles::parse_bytes(b"C/C=C/C=C/C").unwrap().into_table_ir(), 3, Some(StereoCoset::Lit(1)))]
    #[case::hexadiene_ez_first(Smiles::parse_bytes(b"C/C=C/C=C\\C").unwrap().into_table_ir(), 1, Some(StereoCoset::Lit(1)))]
    #[case::hexadiene_ez_second(Smiles::parse_bytes(b"C/C=C/C=C\\C").unwrap().into_table_ir(), 3, Some(StereoCoset::Lit(0)))]
    fn test_raise_cis_trans_stereo(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: Option<StereoCoset>,
    ) {
        let expected = expected
            .map(|coset| BondConstraintForm::CisTransStereo(CisTransStereoForm::stereo(coset)));
        assert_eq!(raise_cis_trans_stereo(&mol, bond_idx), Ok(expected));
    }

    #[rstest]
    #[case::conflict(Smiles::parse_bytes(b"F/C(\\Cl)=CF").unwrap().into_table_ir(), 2, RaiseError::CisTransConflict { atom: 1 })]
    fn test_raise_cis_trans_stereo_error(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: RaiseError,
    ) {
        assert_eq!(raise_cis_trans_stereo(&mol, bond_idx), Err(expected));
    }

    #[rstest]
    #[case::dangling(Smiles::parse_bytes(b"F/C=C").unwrap().into_table_ir(), 0, Err(RaiseError::DanglingBondDirection { bond: 0 }))]
    #[case::flanks_capable(Smiles::parse_bytes(b"O=C1/C=C\\CCC1").unwrap().into_table_ir(), 2, Ok(()))]
    fn test_validate_bond_direction(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: Result<(), RaiseError>,
    ) {
        assert_eq!(validate_bond_direction(&mol, bond_idx), expected);
    }
}
