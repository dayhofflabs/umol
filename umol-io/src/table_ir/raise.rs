//! TableIR → `umol_graph_ir::Molecule` raise.
//!
//! Implements `TryIntoIr<Molecule> for &Molecule` (and the per-atom and
//! per-bond analogues). Table IR fields copy to `Lit` / `Undetermined`; IO
//! raise applies fixed IO ground semantics for resolution.

use std::any::Any;
use std::collections::HashSet;

use thiserror::Error;
use umol_graph_ir::ir::{
    AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId, BondConstraintForm, BondForm,
    BooleanForm, CisTransStereoForm, Constraints, DativeBondForm, ElementForm, IsotopeMassForm,
    Lattice, Molecule, MoleculeEntries, MoleculeIntegrityError, MulticenterBondForm,
    NoncovalentBondForm, NumForm, StereoAtomForm, StereoCoset, StereoKind, StereoLigand,
    StereoLigandKind, TetrahedralStereoForm, TryIntoIr, UnpairedElectronsForm,
};
use umol_perm::{ClassKey, Permutation};
use umol_utils::error::UmolError;

use crate::table_ir::atom::Atom as TableAtom;
use crate::table_ir::bond::{
    Bond as TableBond, BondDonation as TableBondDonation, BondOrder as TableBondOrder,
};
use crate::table_ir::raise::utils::coset_from_wedge_winding;
use crate::table_ir::{
    AtomNeighbors, BondStereo, Chirality, ChiralityFrame, Molecule as TableMolecule,
    StereoLigand as TableStereoLigand, Winding,
};

mod utils;

use utils::{
    cis_trans_capable, cis_trans_side, cis_trans_sides_from_positions, double_bond_partner,
    first_neighbor_toward_ordering, has_either_wedge, neighbor_count, noncovalent_kind,
    tetrahedral_ligand_ordering, validate_bond_direction, validate_tetrahedral_geometry,
    wedge_bond_neighbors, StereoBondAtom, StereoHalfplane,
};

/// Error variants for TableIR -> Molecule raise.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum RaiseError {
    #[error(transparent)]
    MoleculeEntries(#[from] MoleculeIntegrityError),
    #[error("tetrahedral stereo at atom {atom} with {count} ligands, expected 3 or 4 ligands")]
    TetrahedralLigandCount { atom: usize, count: usize },
    #[error("directional bond {bond} not adjacent to a stereogenic double bond")]
    DanglingBondDirection { bond: usize },
    #[error("contradictory cis/trans markers at atom {atom}")]
    CisTransConflict { atom: usize },
    #[error("inconsistent wedge bonds at atom {atom}")]
    WedgeConflict { atom: usize },
    #[error("wedge coordinates at atom {atom} do not determine a configuration")]
    DegenerateWedgeGeometry { atom: usize },
}

impl UmolError for RaiseError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl TryIntoIr<Molecule> for &TableMolecule {
    type Context = ();
    type Error = RaiseError;

    fn try_into_ir(self, context: &Self::Context) -> Result<Molecule, RaiseError> {
        let neighbors = self.atom_neighbors();
        let atoms: Vec<AtomForm> = self
            .atoms
            .iter()
            .enumerate()
            .map(|(atom_idx, table_atom)| {
                let mut atom = table_atom.try_into_ir(context)?;
                if let Some(constraint) = raise_tetrahedral_stereo(self, &neighbors, atom_idx)? {
                    atom.constraints.set(constraint);
                }
                Ok(atom)
            })
            .collect::<Result<_, RaiseError>>()?;

        let mut bonds = Vec::new();
        let mut dative_bonds: Vec<(Vec<AtomId>, AtomId, DativeBondForm)> = Vec::new();
        let mut noncovalent_bonds = Vec::new();
        for (bond_idx, b) in self.bonds.iter().enumerate() {
            validate_bond_direction(self, &neighbors, bond_idx)?;
            let a_idx = AtomId(b.atoms.first());
            let b_idx = AtomId(b.atoms.second());
            if let Some(kind) = b.noncovalent.map(noncovalent_kind) {
                noncovalent_bonds.push(([a_idx, b_idx], NoncovalentBondForm::from_kind(kind)));
            } else if let Some(donation) = b.donation {
                let (donor, acceptor) = match donation {
                    TableBondDonation::Donating => (a_idx, b_idx),
                    TableBondDonation::Accepting => (b_idx, a_idx),
                    _ => {
                        bonds.push((a_idx, b_idx, b.try_into_ir(context)?));
                        continue;
                    }
                };
                let dative_bond = DativeBondForm::new(raise_bond_order(b.order));
                dative_bonds.push((vec![donor], acceptor, dative_bond));
            } else {
                let mut bond_form = b.try_into_ir(context)?;
                if let Some(constraint) = raise_cis_trans_stereo(self, &neighbors, bond_idx)? {
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

        let stereo_atoms = self
            .stereo_atoms
            .iter()
            .map(|frame| {
                let ligands = frame
                    .ligands
                    .iter()
                    .map(|ligand| {
                        let (atom, kind) = match ligand {
                            TableStereoLigand::Atom(atom) => (*atom, StereoLigandKind::Atom),
                            TableStereoLigand::ImplicitHydrogen => {
                                (frame.atom, StereoLigandKind::ImplicitHydrogen)
                            }
                            TableStereoLigand::LonePair => (frame.atom, StereoLigandKind::LonePair),
                        };
                        StereoLigand::new(AtomId(atom), kind)
                    })
                    .collect();
                let coset = match frame.winding {
                    Winding::CounterClockwise => 0,
                    Winding::Clockwise => 1,
                };
                (
                    AtomId(frame.atom),
                    ligands,
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(coset)),
                )
            })
            .collect();
        let constraints = Constraints::new();

        Molecule::try_from_entries(MoleculeEntries {
            atoms,
            bonds,
            stereo_atoms,
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
    type Context = ();
    type Error = RaiseError;

    fn try_into_ir(self, _context: &Self::Context) -> Result<AtomForm, RaiseError> {
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
        // Closed shell is implied only where the notation owns the hydrogen
        // count (implicit-valence atoms); an explicit hydrogen count leaves
        // the shell open for resolution ([H], [CH3] resolve as radicals).
        if matches!(atom.unpaired_electrons.count, NumForm::Undetermined)
            && self.implicit_hydrogens.is_none()
        {
            atom.unpaired_electrons.count = NumForm::Lit(0);
        }
        atom.constraints.retain(|c| !c.is_undetermined());
        Ok(atom)
    }
}

impl TryIntoIr<BondForm> for &TableBond {
    type Context = ();
    type Error = RaiseError;

    fn try_into_ir(self, _context: &Self::Context) -> Result<BondForm, RaiseError> {
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
    neighbors: &AtomNeighbors,
    atom_idx: usize,
) -> Result<Option<AtomConstraintForm>, RaiseError> {
    if mol
        .stereo_atoms
        .iter()
        .any(|frame| frame.atom as usize == atom_idx)
    {
        return Ok(None);
    }
    // CTfile atom parity is retained in TableIR and not read: the specification marks the field
    // ignored when read, and a CTfile record's tetrahedral stereo comes from its wedges.
    let chirality = match mol.chirality_frame {
        Some(ChiralityFrame::FirstNeighborToward) => mol.atoms[atom_idx].chirality,
        Some(ChiralityFrame::LastNeighborAway) | None => None,
    };
    let (relabeling, source_coset): (Permutation, usize) = match chirality {
        Some(symbol) => {
            let source_coset = match symbol {
                Chirality::CounterClockwise | Chirality::Tetrahedral { arr: 1 } => 0,
                Chirality::Clockwise | Chirality::Tetrahedral { arr: 2 } => 1,
                _ => return Ok(None),
            };
            validate_tetrahedral_geometry(neighbors, atom_idx)?;
            let permutation = Permutation::between(
                &first_neighbor_toward_ordering(neighbors, atom_idx),
                &tetrahedral_ligand_ordering(neighbors, atom_idx),
            )
            .expect("validated tetrahedral frames contain the same ligands");
            (permutation, source_coset)
        }
        None => {
            let wedged = wedge_bond_neighbors(mol, neighbors, atom_idx);
            // An either wedge at an atom with one double bond marks that bond; see
            // `raise_cis_trans_stereo`.
            if has_either_wedge(mol, neighbors, atom_idx)
                && double_bond_partner(mol, neighbors, atom_idx).is_none()
            {
                if !wedged.is_empty() {
                    return Err(RaiseError::WedgeConflict { atom: atom_idx });
                }
                validate_tetrahedral_geometry(neighbors, atom_idx)?;
                return Ok(Some(AtomConstraintForm::TetrahedralStereo(
                    TetrahedralStereoForm::stereo(StereoCoset::Undetermined),
                )));
            }
            let Some(positions) = mol.positions.as_ref() else {
                return Ok(None);
            };
            let count = neighbor_count(neighbors, atom_idx);
            if count != 3 && count != 4 {
                return Ok(None);
            }
            let Some(&(neighbor_idx, outofplane)) = wedged.first() else {
                return Ok(None);
            };
            let target_ordering = tetrahedral_ligand_ordering(neighbors, atom_idx);
            let winding = |neighbor_idx, outofplane| {
                coset_from_wedge_winding(
                    atom_idx,
                    &target_ordering,
                    neighbor_idx,
                    positions,
                    outofplane,
                )
                .ok_or(RaiseError::DegenerateWedgeGeometry { atom: atom_idx })
            };
            let source_coset = winding(neighbor_idx, outofplane)?;
            for &(neighbor_idx, outofplane) in &wedged[1..] {
                if winding(neighbor_idx, outofplane)? != source_coset {
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
    neighbors: &AtomNeighbors,
    bond_idx: usize,
) -> Result<Option<BondConstraintForm>, RaiseError> {
    let bond = &mol.bonds[bond_idx];
    if bond.order != TableBondOrder::Double {
        return Ok(None);
    }
    let atom_1_idx = bond.start_atom() as usize;
    let atom_2_idx = bond.end_atom() as usize;
    // Stereo code 3, or the drawing convention of an either wedge at an atom of the double bond,
    // asserts an unknown configuration.
    let either_marked = |atom_idx: usize, other_atom_idx: usize| {
        has_either_wedge(mol, neighbors, atom_idx)
            && double_bond_partner(mol, neighbors, atom_idx) == Some(other_atom_idx)
    };
    if bond.stereo == Some(BondStereo::Either)
        || either_marked(atom_1_idx, atom_2_idx)
        || either_marked(atom_2_idx, atom_1_idx)
    {
        return Ok(Some(BondConstraintForm::CisTransStereo(
            CisTransStereoForm::stereo(StereoCoset::Undetermined),
        )));
    }
    if !cis_trans_capable(neighbors, atom_1_idx, atom_2_idx) {
        return Ok(None);
    }
    // Directional marks decide when present at both atoms; a double bond without any mark is read
    // from coordinates when the record has them.
    let (side_1, side_2) = match (
        cis_trans_side(mol, neighbors, atom_1_idx, atom_2_idx)?,
        cis_trans_side(mol, neighbors, atom_2_idx, atom_1_idx)?,
    ) {
        (Some(side_1), Some(side_2)) => (side_1, side_2),
        (None, None) => {
            let Some(positions) = mol.positions.as_ref() else {
                return Ok(None);
            };
            match cis_trans_sides_from_positions(mol, neighbors, atom_1_idx, atom_2_idx, positions)
            {
                Some((side_1, side_2)) => (side_1, side_2),
                None => return Ok(None),
            }
        }
        _ => return Ok(None),
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
    use umol_graph_ir::ir::{AtomConstraintsForm, BondId, Entity, StereoAtomId};

    use super::*;
    use crate::ctfile::parse_mol_to_ir;
    use crate::ctfile::parser::parse_mol_bytes_to_table_ir;
    use crate::smiles::Smiles;
    use crate::smiles::SmilesIoConfig;
    use crate::table_ir::atom::Atom as TableAtom;
    use crate::table_ir::bond::{Bond as TableBond, BondOrder as TableBondOrder};
    use crate::table_ir::{Molecule as TableMolecule, StereoAtom};

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

    const CFCLBRI_EITHER_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.6906   -0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.6906    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.6906    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  0        0\n  2  4  1  0        0\n  2  5  1  4        0\nM  END\n";

    const CFCLBRI_MIXED_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.6906   -0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.6906    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.6906    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  1        0\n  2  4  1  0        0\n  2  5  1  4        0\nM  END\n";

    const EITHER_WEDGE_TWO_LIGANDS_MOL: &str = "\n\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n   -0.6906    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6906    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  4        0\n  2  3  1  0        0\nM  END\n";

    // Atom 1 is wedged toward atom 4, whose four neighbors make it look like a stereo site.
    const WIDE_ENDPOINT_WEDGE_MOL: &str = "\n\n\n  8  7  0  0  1  0  0  0  0  0999 V2000\n    0.6906   -0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.6906    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.6906    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -1.3812    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6906   -1.3812    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.3812    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  0        0\n  2  4  1  0        0\n  2  5  1  1        0\n  5  6  1  0        0\n  5  7  1  0        0\n  5  8  1  0        0\nM  END\n";

    // Atoms 1 and 2 are each wedged toward atom 0, written as the higher-numbered first atom.
    const SHARED_WIDE_ENDPOINT_WEDGE_MOL: &str = "\n\n\n 10  9  0  0  1  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000   -1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n   -2.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0000    1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0000   -1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1        0\n  3  1  1  6        0\n  1  4  1  0        0\n  2  5  1  0        0\n  2  6  1  0        0\n  2  7  1  0        0\n  3  8  1  0        0\n  3  9  1  0        0\n  3 10  1  0        0\nM  END\n";

    // Atom 0 is wedged toward atom 1, which has its own wedge toward atom 4; read alone, the
    // incoming wedge would give atom 1 the opposite coset.
    const INCOMING_WEDGE_MOL: &str = "\n\n\n  8  7  0  0  1  0  0  0  0  0999 V2000\n    0.6906   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.6906    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906   -0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.6906    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n    1.3812   -0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6906    0.6906    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.6906   -0.6906    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  1        0\n  2  3  1  0        0\n  2  4  1  0        0\n  2  5  1  1        0\n  1  6  1  0        0\n  1  7  1  0        0\n  1  8  1  0        0\nM  END\n";

    const DIFLUOROETHENE_E_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500   -0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  2  0        0\n  3  4  1  0        0\nM  END\n";
    const DIFLUOROETHENE_E_REVERSED_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500   -0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  3  2  2  0        0\n  3  4  1  0        0\nM  END\n";

    const DIFLUOROETHENE_Z_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500    0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  2  0        0\n  3  4  1  0        0\nM  END\n";

    const DIFLUOROETHENE_E_3D_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.0000    0.7500 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500    0.0000   -0.7500 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  2  0        0\n  3  4  1  0        0\nM  END\n";

    const DIFLUOROETHENE_Z_3D_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.0000    0.7500 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500    0.0000    0.7500 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  2  0        0\n  3  4  1  0        0\nM  END\n";

    const CYCLOHEXENE_MOL: &str = "\n\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5000    0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5000   -0.8660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  6  1  0        0\n  1  2  2  0        0\n  2  3  1  0        0\n  3  4  1  0        0\n  4  5  1  0        0\n  5  6  1  0        0\nM  END\n";

    const DIFLUOROETHENE_ZERO_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  2  0        0\n  3  4  1  0        0\nM  END\n";

    const DIFLUOROETHENE_COLLINEAR_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500   -0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  2  0        0\n  3  4  1  0        0\nM  END\n";
    const DIFLUOROETHENE_COLLINEAR_REVERSED_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500   -0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  3  2  2  0        0\n  3  4  1  0        0\nM  END\n";

    const CFCLBRI_ZERO_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  0        0\n  2  4  1  0        0\n  2  5  1  1        0\nM  END\n";

    const CFCLBRI_COLLINEAR_WEDGE_MOL: &str = "\n\n\n  5  4  0  0  1  0  0  0  0  0999 V2000\n    0.6906    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3812    0.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6906    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3812    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  1  0        0\n  2  4  1  0        0\n  2  5  1  1        0\nM  END\n";

    const FOLDED_ALKENE_MOL: &str = "\n\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6500    0.2500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500   -0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0        0\n  2  3  1  0        0\n  3  4  2  0        0\n  4  5  1  0        0\nM  END\n";

    const WAVY_ALKENE_TWO_LIGANDS_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  4        0\n  2  3  2  0        0\n  3  4  1  0        0\nM  END\n";

    const WAVY_ALKENE_THREE_LIGANDS_MOL: &str = "\n\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9500   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6500   -0.7500    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  4        0\n  2  3  2  0        0\n  3  4  1  0        0\n  2  5  1  0        0\nM  END\n";

    const ALLENE_MOL: &str = "\n\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2500   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  2  0        0\n  3  4  2  0        0\n  4  5  1  0        0\nM  END\n";

    const ISOTHIOCYANATE_MOL: &str = "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -0.6500    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.3000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.6000    0.0000    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0        0\n  2  3  2  0        0\n  3  4  2  0        0\nM  END\n";

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
                        count: NumForm::Undetermined,
                        multiplicity: NumForm::Undetermined,
                    },
                    constraints: AtomConstraintsForm::new(),
                }],
                ..Default::default()
            })
        );
    }

    #[rstest]
    #[case::actual(
        Smiles::parse("[C@](F)(Cl)(Br)I").unwrap().into_table_ir(), 0,
        vec![TableStereoLigand::Atom(4), TableStereoLigand::Atom(2), TableStereoLigand::Atom(1), TableStereoLigand::Atom(3)], Winding::Clockwise,
        vec![(4, StereoLigandKind::Atom), (2, StereoLigandKind::Atom), (1, StereoLigandKind::Atom), (3, StereoLigandKind::Atom)], 1
    )]
    #[case::hydrogen(
        Smiles::parse("C.[C@@H](F)(Cl)Br").unwrap().into_table_ir(), 1,
        vec![TableStereoLigand::Atom(4), TableStereoLigand::ImplicitHydrogen, TableStereoLigand::Atom(2), TableStereoLigand::Atom(3)], Winding::CounterClockwise,
        vec![(4, StereoLigandKind::Atom), (1, StereoLigandKind::ImplicitHydrogen), (2, StereoLigandKind::Atom), (3, StereoLigandKind::Atom)], 0
    )]
    #[case::lone_pair(
        Smiles::parse("[N@@](C)(F)Cl").unwrap().into_table_ir(), 0,
        vec![TableStereoLigand::Atom(1), TableStereoLigand::LonePair, TableStereoLigand::Atom(3), TableStereoLigand::Atom(2)], Winding::CounterClockwise,
        vec![(1, StereoLigandKind::Atom), (0, StereoLigandKind::LonePair), (3, StereoLigandKind::Atom), (2, StereoLigandKind::Atom)], 0
    )]
    #[case::mixed_wedges(
        parse_mol_bytes_to_table_ir(CFCLBRI_MIXED_WEDGE_MOL.as_bytes()).unwrap(), 1,
        vec![TableStereoLigand::Atom(4), TableStereoLigand::Atom(0), TableStereoLigand::Atom(2), TableStereoLigand::Atom(3)], Winding::Clockwise,
        vec![(4, StereoLigandKind::Atom), (0, StereoLigandKind::Atom), (2, StereoLigandKind::Atom), (3, StereoLigandKind::Atom)], 1
    )]
    fn test_table_molecule_try_into_ir_frame(
        #[case] mut table: TableMolecule,
        #[case] atom: u32,
        #[case] ligands: Vec<TableStereoLigand>,
        #[case] winding: Winding,
        #[case] expected_ligands: Vec<(u32, StereoLigandKind)>,
        #[case] coset: u32,
    ) {
        table.stereo_atoms = vec![StereoAtom {
            atom,
            ligands,
            winding,
        }];
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: table
                .atoms
                .iter()
                .map(|atom| atom.try_into_ir(&()).unwrap())
                .collect(),
            bonds: table
                .bonds
                .iter()
                .map(|bond| {
                    (
                        AtomId(bond.start_atom()),
                        AtomId(bond.end_atom()),
                        bond.try_into_ir(&()).unwrap(),
                    )
                })
                .collect(),
            stereo_atoms: vec![(
                AtomId(atom),
                expected_ligands
                    .into_iter()
                    .map(|(atom, kind)| StereoLigand::new(AtomId(atom), kind))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(coset)),
            )],
            ..Default::default()
        });
        assert_eq!((&table).try_into_ir(&()), Ok(expected));
    }

    #[rstest]
    #[case::empty(1, vec![], 1, MoleculeIntegrityError::StereoLigandArity { entity: Entity::StereoAtom(StereoAtomId(0)), kind: StereoKind::Tetrahedral, expected: 4, actual: 0 })]
    #[case::oversized(1, vec![TableStereoLigand::Atom(0), TableStereoLigand::Atom(2), TableStereoLigand::Atom(3), TableStereoLigand::Atom(4), TableStereoLigand::ImplicitHydrogen], 1, MoleculeIntegrityError::StereoLigandArity { entity: Entity::StereoAtom(StereoAtomId(0)), kind: StereoKind::Tetrahedral, expected: 4, actual: 5 })]
    #[case::repeated_hydrogen(1, vec![TableStereoLigand::ImplicitHydrogen, TableStereoLigand::ImplicitHydrogen, TableStereoLigand::Atom(2), TableStereoLigand::Atom(3)], 1, MoleculeIntegrityError::DuplicateStereoLigand { entity: Entity::StereoAtom(StereoAtomId(0)), ligand: StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen) })]
    #[case::missing_site(8, vec![TableStereoLigand::Atom(0), TableStereoLigand::Atom(2), TableStereoLigand::Atom(3), TableStereoLigand::Atom(4)], 1, MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(8)) })]
    #[case::missing_ligand(1, vec![TableStereoLigand::Atom(8), TableStereoLigand::Atom(2), TableStereoLigand::Atom(3), TableStereoLigand::Atom(4)], 1, MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(8)) })]
    #[case::duplicate_site(1, vec![TableStereoLigand::Atom(0), TableStereoLigand::Atom(2), TableStereoLigand::Atom(3), TableStereoLigand::Atom(4)], 2, MoleculeIntegrityError::StereoAtomSitesDuplicate { atom: AtomId(1) })]
    #[case::site_ligand(1, vec![TableStereoLigand::Atom(1), TableStereoLigand::Atom(2), TableStereoLigand::Atom(3), TableStereoLigand::Atom(4)], 1, MoleculeIntegrityError::DuplicateParticipant { entity: Entity::StereoAtom(StereoAtomId(0)), atom: AtomId(1) })]
    #[case::non_neighbor(1, vec![TableStereoLigand::Atom(5), TableStereoLigand::Atom(2), TableStereoLigand::Atom(3), TableStereoLigand::Atom(4)], 1, MoleculeIntegrityError::StereoLigandIncidenceMismatch { entity: Entity::StereoAtom(StereoAtomId(0)) })]
    fn test_table_molecule_try_into_ir_frame_error(
        #[case] atom: u32,
        #[case] ligands: Vec<TableStereoLigand>,
        #[case] copies: usize,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut table = parse_mol_bytes_to_table_ir(CFCLBRI_MIXED_WEDGE_MOL.as_bytes()).unwrap();
        table.atoms.push(TableAtom::from_element(Element::C));
        table.stereo_atoms = vec![
            StereoAtom {
                atom,
                ligands,
                winding: Winding::Clockwise
            };
            copies
        ];
        if atom != 1 {
            table.bonds.iter_mut().for_each(|bond| bond.wedge = None);
        }
        let result: Result<Molecule, _> = (&table).try_into_ir(&());
        assert_eq!(result, Err(RaiseError::MoleculeEntries(expected)));
    }

    #[rstest]
    fn test_table_molecule_try_into_ir_parity() {
        let mut table = parse_mol_bytes_to_table_ir(METHANE_MOL.as_bytes()).unwrap();
        let expected: Molecule = (&table).try_into_ir(&()).unwrap();
        table.atoms[0].chirality = Some(Chirality::Clockwise);
        table.chirality_frame = Some(ChiralityFrame::LastNeighborAway);
        assert_eq!((&table).try_into_ir(&()), Ok(expected));
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
        RaiseError::MoleculeEntries(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        })
    )]
    #[case::repeated_virtual_tetrahedral_completion(
        Smiles::parse_bytes(b"[C@H2](F)Cl").unwrap().into_table_ir(),
        RaiseError::MoleculeEntries(MoleculeIntegrityError::DuplicateStereoLigand {
            entity: Entity::StereoAtom(StereoAtomId(0)),
            ligand: StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        })
    )]
    #[case::incomplete(
        Smiles::parse("[C@]").unwrap().into_table_ir(),
        RaiseError::MoleculeEntries(MoleculeIntegrityError::StereoLigandArity { entity: Entity::StereoAtom(StereoAtomId(0)), kind: StereoKind::Tetrahedral, expected: 4, actual: 0 })
    )]
    #[case::oversized(
        Smiles::parse("[C@](F)(Cl)(Br)(I)N").unwrap().into_table_ir(),
        RaiseError::MoleculeEntries(MoleculeIntegrityError::StereoLigandArity { entity: Entity::StereoAtom(StereoAtomId(0)), kind: StereoKind::Tetrahedral, expected: 4, actual: 5 })
    )]
    fn test_table_molecule_try_into_ir_error(
        #[case] molecule: TableMolecule,
        #[case] expected: RaiseError,
    ) {
        let actual: Result<Molecule, RaiseError> = (&molecule).try_into_ir(&());
        assert_eq!(actual, Err(expected));
    }

    #[rstest]
    #[case::tetrahedral("C[C@H](N)O", Entity::Atom(AtomId(1)))]
    #[case::cis_trans("F/C=C/F", Entity::Bond(BondId(1)))]
    fn test_table_molecule_try_into_ir_stereo(#[case] input: &str, #[case] entity: Entity) {
        let smiles = Smiles::parse(input).unwrap();
        let molecule: Molecule = smiles.as_table_ir().try_into_ir(&()).unwrap();

        assert_eq!(molecule.stereo_bonds().count(), 0);
        match entity {
            Entity::Atom(id) => {
                let frames: Vec<_> = molecule
                    .stereo_atoms()
                    .iter()
                    .map(|frame| {
                        (
                            frame.site_id(),
                            frame
                                .ligands()
                                .map(|ligand| (ligand.atom_id(), ligand.kind()))
                                .collect::<Vec<_>>(),
                            frame.attributes.clone(),
                        )
                    })
                    .collect();
                assert_eq!(
                    frames,
                    vec![(
                        id,
                        vec![
                            (AtomId(0), StereoLigandKind::Atom),
                            (id, StereoLigandKind::ImplicitHydrogen),
                            (AtomId(2), StereoLigandKind::Atom),
                            (AtomId(3), StereoLigandKind::Atom)
                        ],
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                    )]
                );
                assert_eq!(
                    molecule
                        .atom(id)
                        .attributes
                        .constraints
                        .tetrahedral_stereo(),
                    None
                );
            }
            Entity::Bond(id) => assert!(molecule
                .bond(id)
                .attributes
                .constraints
                .cis_trans_stereo()
                .is_some()),
            _ => unreachable!("stereo constraints are raised only on atoms and bonds"),
        }
    }

    #[rstest]
    #[case::tetrahedral(CFCLBRI_SINGLE_WEDGE_MOL, Entity::Atom(AtomId(1)))]
    #[case::cis_trans(CIS_TRANS_EITHER_MOL, Entity::Bond(BondId(1)))]
    fn test_parse_mol_to_ir_stereo(#[case] input: &str, #[case] entity: Entity) {
        let molecule = parse_mol_to_ir(input).unwrap();

        assert_eq!(molecule.stereo_atoms().count(), 0);
        assert_eq!(molecule.stereo_bonds().count(), 0);
        match entity {
            Entity::Atom(id) => assert!(molecule
                .atom(id)
                .attributes
                .constraints
                .tetrahedral_stereo()
                .is_some()),
            Entity::Bond(id) => assert!(molecule
                .bond(id)
                .attributes
                .constraints
                .cis_trans_stereo()
                .is_some()),
            _ => unreachable!("stereo constraints are raised only on atoms and bonds"),
        }
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

    // A bare aromatic atom stays H-open — the hydrogen count is the valence
    // model's decision; explicit bracket H is preserved.
    #[rstest]
    #[case::aromatic_nitrogen_bare(Element::N, Some(true), None, NumForm::Undetermined)]
    #[case::aromatic_oxygen_bare(Element::O, Some(true), None, NumForm::Undetermined)]
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
    #[case::carbon_h0(CARBON_H0_EXPLICIT_MOL, "C#i=#c0#h0")]
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
    #[case::mol_wedge_cfclbri(parse_mol_bytes_to_table_ir(CFCLBRI_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_cfclbri_r(parse_mol_bytes_to_table_ir(CFCLBRI_R_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_wedge_cfclbri_single(parse_mol_bytes_to_table_ir(CFCLBRI_SINGLE_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_wedge_butanol(parse_mol_bytes_to_table_ir(BUTANOL_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_alanine(parse_mol_bytes_to_table_ir(ALANINE_WEDGE_MOL.as_bytes()).unwrap(), 4, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_sulfoxide(parse_mol_bytes_to_table_ir(SULFOXIDE_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_methyloxirane(parse_mol_bytes_to_table_ir(METHYLOXIRANE_WEDGE_MOL.as_bytes()).unwrap(), 0, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_prochiral_methylene(parse_mol_bytes_to_table_ir(PROCHIRAL_METHYLENE_WEDGE_MOL.as_bytes()).unwrap(), 0, Some(StereoCoset::Lit(1)))]
    #[case::mol_wedge_wide_endpoint_site(parse_mol_bytes_to_table_ir(WIDE_ENDPOINT_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_wedge_wide_endpoint(parse_mol_bytes_to_table_ir(WIDE_ENDPOINT_WEDGE_MOL.as_bytes()).unwrap(), 4, None)]
    #[case::mol_wedge_shared_wide_endpoint_first_site(parse_mol_bytes_to_table_ir(SHARED_WIDE_ENDPOINT_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_shared_wide_endpoint_second_site(parse_mol_bytes_to_table_ir(SHARED_WIDE_ENDPOINT_WEDGE_MOL.as_bytes()).unwrap(), 2, Some(StereoCoset::Lit(0)))]
    #[case::mol_wedge_shared_wide_endpoint(parse_mol_bytes_to_table_ir(SHARED_WIDE_ENDPOINT_WEDGE_MOL.as_bytes()).unwrap(), 0, None)]
    #[case::mol_wedge_incoming_site(parse_mol_bytes_to_table_ir(INCOMING_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_wedge_incoming_source(parse_mol_bytes_to_table_ir(INCOMING_WEDGE_MOL.as_bytes()).unwrap(), 0, Some(StereoCoset::Lit(0)))]
    #[case::mol_either_wedge(parse_mol_bytes_to_table_ir(CFCLBRI_EITHER_WEDGE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Undetermined))]
    #[case::mol_either_wedge_wide_endpoint(parse_mol_bytes_to_table_ir(CFCLBRI_EITHER_WEDGE_MOL.as_bytes()).unwrap(), 4, None)]
    #[case::cx_wiggly(Smiles::parse_bytes_with(b"F[C](Cl)(Br)I |w:1.0|", &SmilesIoConfig::chemaxon()).unwrap().into_table_ir(), 1, Some(StereoCoset::Undetermined))]
    #[case::cx_wiggly_up(Smiles::parse_bytes_with(b"F[C](Cl)(Br)I |wU:1.0|", &SmilesIoConfig::chemaxon()).unwrap().into_table_ir(), 1, Some(StereoCoset::Undetermined))]
    #[case::cx_wiggly_down(Smiles::parse_bytes_with(b"F[C](Cl)(Br)I |wD:1.0|", &SmilesIoConfig::chemaxon()).unwrap().into_table_ir(), 1, Some(StereoCoset::Undetermined))]
    #[case::cx_wiggly_wide_endpoint(Smiles::parse_bytes_with(b"F[C](Cl)(Br)I |w:1.0|", &SmilesIoConfig::chemaxon()).unwrap().into_table_ir(), 0, None)]
    #[case::mol_wavy_alkene_two_ligands(parse_mol_bytes_to_table_ir(WAVY_ALKENE_TWO_LIGANDS_MOL.as_bytes()).unwrap(), 1, None)]
    #[case::mol_wavy_alkene_three_ligands(parse_mol_bytes_to_table_ir(WAVY_ALKENE_THREE_LIGANDS_MOL.as_bytes()).unwrap(), 1, None)]
    #[case::cx_wiggly_alkene(Smiles::parse_bytes_with(b"CC=CC |w:1.0|", &SmilesIoConfig::chemaxon()).unwrap().into_table_ir(), 1, None)]
    #[case::no_descriptor(Smiles::parse_bytes(b"F[C@](Cl)(Br)I").unwrap().into_table_ir(), 0, None)]
    fn test_raise_tetrahedral_stereo(
        #[case] mol: TableMolecule,
        #[case] atom_idx: usize,
        #[case] expected: Option<StereoCoset>,
    ) {
        let expected = expected.map(|coset| {
            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::stereo(coset))
        });
        assert_eq!(
            raise_tetrahedral_stereo(&mol, &mol.atom_neighbors(), atom_idx),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::wedge_conflict(parse_mol_bytes_to_table_ir(WEDGE_CONFLICT_MOL.as_bytes()).unwrap(), 0, RaiseError::WedgeConflict { atom: 0 })]
    #[case::cfclbri_inconsistent_wedges(parse_mol_bytes_to_table_ir(CFCLBRI_INCONSISTENT_WEDGE_MOL.as_bytes()).unwrap(), 1, RaiseError::WedgeConflict { atom: 1 })]
    #[case::cfclbri_definite_and_either_wedge(parse_mol_bytes_to_table_ir(CFCLBRI_MIXED_WEDGE_MOL.as_bytes()).unwrap(), 1, RaiseError::WedgeConflict { atom: 1 })]
    #[case::either_wedge_two_ligands(parse_mol_bytes_to_table_ir(EITHER_WEDGE_TWO_LIGANDS_MOL.as_bytes()).unwrap(), 1, RaiseError::TetrahedralLigandCount { atom: 1, count: 2 })]
    #[case::wedge_zero_coordinates(parse_mol_bytes_to_table_ir(CFCLBRI_ZERO_WEDGE_MOL.as_bytes()).unwrap(), 1, RaiseError::DegenerateWedgeGeometry { atom: 1 })]
    #[case::wedge_collinear_coordinates(parse_mol_bytes_to_table_ir(CFCLBRI_COLLINEAR_WEDGE_MOL.as_bytes()).unwrap(), 1, RaiseError::DegenerateWedgeGeometry { atom: 1 })]
    fn test_raise_tetrahedral_stereo_error(
        #[case] mol: TableMolecule,
        #[case] atom_idx: usize,
        #[case] expected: RaiseError,
    ) {
        assert_eq!(
            raise_tetrahedral_stereo(&mol, &mol.atom_neighbors(), atom_idx),
            Err(expected)
        );
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
    #[case::mol_coordinates_e(parse_mol_bytes_to_table_ir(DIFLUOROETHENE_E_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_coordinates_e_reversed_bond_line(parse_mol_bytes_to_table_ir(DIFLUOROETHENE_E_REVERSED_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_coordinates_z(parse_mol_bytes_to_table_ir(DIFLUOROETHENE_Z_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_coordinates_e_3d(parse_mol_bytes_to_table_ir(DIFLUOROETHENE_E_3D_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(1)))]
    #[case::mol_coordinates_z_3d(parse_mol_bytes_to_table_ir(DIFLUOROETHENE_Z_3D_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_coordinates_ring(parse_mol_bytes_to_table_ir(CYCLOHEXENE_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Lit(0)))]
    #[case::mol_coordinates_folded(parse_mol_bytes_to_table_ir(FOLDED_ALKENE_MOL.as_bytes()).unwrap(), 2, None)]
    #[case::mol_coordinates_zero(parse_mol_bytes_to_table_ir(DIFLUOROETHENE_ZERO_MOL.as_bytes()).unwrap(), 1, None)]
    #[case::mol_coordinates_collinear_substituent(parse_mol_bytes_to_table_ir(DIFLUOROETHENE_COLLINEAR_MOL.as_bytes()).unwrap(), 1, None)]
    #[case::mol_coordinates_collinear_substituent_reversed_bond_line(parse_mol_bytes_to_table_ir(DIFLUOROETHENE_COLLINEAR_REVERSED_MOL.as_bytes()).unwrap(), 1, None)]
    #[case::mol_wavy_two_ligands(parse_mol_bytes_to_table_ir(WAVY_ALKENE_TWO_LIGANDS_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Undetermined))]
    #[case::mol_wavy_three_ligands(parse_mol_bytes_to_table_ir(WAVY_ALKENE_THREE_LIGANDS_MOL.as_bytes()).unwrap(), 1, Some(StereoCoset::Undetermined))]
    #[case::cx_wiggly_alkene(Smiles::parse_bytes_with(b"CC=CC |w:1.0|", &SmilesIoConfig::chemaxon()).unwrap().into_table_ir(), 1, Some(StereoCoset::Undetermined))]
    #[case::mol_allene_first(parse_mol_bytes_to_table_ir(ALLENE_MOL.as_bytes()).unwrap(), 1, None)]
    #[case::mol_allene_second(parse_mol_bytes_to_table_ir(ALLENE_MOL.as_bytes()).unwrap(), 2, None)]
    #[case::mol_isothiocyanate_n_c(parse_mol_bytes_to_table_ir(ISOTHIOCYANATE_MOL.as_bytes()).unwrap(), 1, None)]
    #[case::mol_isothiocyanate_c_s(parse_mol_bytes_to_table_ir(ISOTHIOCYANATE_MOL.as_bytes()).unwrap(), 2, None)]
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
        assert_eq!(
            raise_cis_trans_stereo(&mol, &mol.atom_neighbors(), bond_idx),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::conflict(Smiles::parse_bytes(b"F/C(\\Cl)=CF").unwrap().into_table_ir(), 2, RaiseError::CisTransConflict { atom: 1 })]
    fn test_raise_cis_trans_stereo_error(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: RaiseError,
    ) {
        assert_eq!(
            raise_cis_trans_stereo(&mol, &mol.atom_neighbors(), bond_idx),
            Err(expected)
        );
    }

    #[rstest]
    #[case::dangling(Smiles::parse_bytes(b"F/C=C").unwrap().into_table_ir(), 0, Err(RaiseError::DanglingBondDirection { bond: 0 }))]
    #[case::flanks_capable(Smiles::parse_bytes(b"O=C1/C=C\\CCC1").unwrap().into_table_ir(), 2, Ok(()))]
    fn test_validate_bond_direction(
        #[case] mol: TableMolecule,
        #[case] bond_idx: usize,
        #[case] expected: Result<(), RaiseError>,
    ) {
        assert_eq!(
            validate_bond_direction(&mol, &mol.atom_neighbors(), bond_idx),
            expected
        );
    }
}
