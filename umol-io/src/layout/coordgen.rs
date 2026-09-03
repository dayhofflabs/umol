//! Graph-IR projection for the CoordGen molecule-layout backend.

use umol_coordgen_sys::{generate_coordinates, Bond, CisTransBond, CoordgenError, SideRelation};
use umol_graph_ir::ir::{
    AsLit, AtomId, Molecule, StereoBondView, StereoCoset, StereoKind, StereoLigand,
    StereoLigandKind,
};

use super::{MoleculeLayout, Point2D};

// CoordGen's fixed working scale is defined by BONDLENGTH in sketcherMinimizerMaths.h.
const COORDGEN_BOND_LENGTH: f64 = 50.0;

pub(crate) fn layout(molecule: &Molecule) -> Result<MoleculeLayout, CoordgenError> {
    let atomic_numbers = molecule
        .atoms()
        .iter()
        .map(|atom| {
            atom.element()
                .as_lit()
                .map(|element| u16::from(element.atomic_number()))
                .unwrap_or(0)
        })
        .collect::<Vec<_>>();
    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            let [atom_0, atom_1] = bond.atom_ids();
            Bond {
                atom_0: atom_0.index(),
                atom_1: atom_1.index(),
                order: match bond.order().as_lit() {
                    Some(order @ 1..=3) => order as u8,
                    _ => 1,
                },
            }
        })
        .collect::<Vec<_>>();
    let cis_trans_bonds = molecule
        .stereo_bonds()
        .iter()
        .filter_map(|stereo| cis_trans_bond(molecule, stereo))
        .collect::<Vec<_>>();

    let positions = generate_coordinates(&atomic_numbers, &bonds, &cis_trans_bonds)?
        .into_iter()
        .map(|point| {
            Point2D::new(
                point.x / COORDGEN_BOND_LENGTH,
                point.y / COORDGEN_BOND_LENGTH,
            )
        })
        .collect();

    Ok(MoleculeLayout::try_new(positions)
        .expect("scaling finite CoordGen points by its finite bond length preserves finiteness"))
}

fn cis_trans_bond(molecule: &Molecule, stereo: StereoBondView<'_>) -> Option<CisTransBond> {
    if stereo.attributes.configuration.kind() != Some(StereoKind::CisTrans) {
        return None;
    }
    let ligands = stereo
        .ligands()
        .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
        .collect::<Vec<_>>();
    let [first_0, first_1, second_0, second_1] = ligands.as_slice() else {
        return None;
    };
    let [site_0, site_1] = stereo.site().atom_ids();
    let (first_pair, second_pair) = if [first_0, first_1]
        .into_iter()
        .all(|ligand| ligand_matches_endpoint(molecule, *ligand, site_0, site_1))
    {
        ([*first_0, *first_1], [*second_0, *second_1])
    } else {
        ([*second_0, *second_1], [*first_0, *first_1])
    };
    let (first_ligand, first_other) = select_actual_ligand(first_pair)?;
    let (second_ligand, second_other) = select_actual_ligand(second_pair)?;
    let requested = [first_ligand, first_other, second_ligand, second_other];
    let StereoCoset::Lit(coset) = stereo.coset_for(requested)? else {
        return None;
    };
    let relation = match coset {
        0 => SideRelation::SameSide,
        1 => SideRelation::OppositeSide,
        _ => return None,
    };

    Some(CisTransBond {
        bond: stereo.site_id().index(),
        first_ligand: first_ligand.atom_id.index(),
        second_ligand: second_ligand.atom_id.index(),
        relation,
    })
}

fn ligand_matches_endpoint(
    molecule: &Molecule,
    ligand: StereoLigand,
    endpoint: AtomId,
    other_endpoint: AtomId,
) -> bool {
    match ligand.kind {
        StereoLigandKind::Atom => {
            ligand.atom_id != other_endpoint
                && molecule
                    .neighbors(endpoint)
                    .any(|neighbor| neighbor.atom_id() == ligand.atom_id)
        }
        StereoLigandKind::ImplicitHydrogen | StereoLigandKind::LonePair => {
            ligand.atom_id == endpoint
        }
    }
}

fn select_actual_ligand(pair: [StereoLigand; 2]) -> Option<(StereoLigand, StereoLigand)> {
    match pair {
        [first, second] if first.kind == StereoLigandKind::Atom => Some((first, second)),
        [first, second] if second.kind == StereoLigandKind::Atom => Some((second, first)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{
        AtomForm, BondForm, BondId, MoleculeEntries, StereoBondForm, StereoTerm,
    };

    use super::*;

    #[rstest]
    #[case::z(StereoCoset::Lit(0), SideRelation::SameSide)]
    #[case::e(StereoCoset::Lit(1), SideRelation::OppositeSide)]
    fn test_cis_trans_bond(#[case] coset: StereoCoset, #[case] relation: SideRelation) {
        let molecule = stereo_molecule(
            StereoKind::CisTrans,
            coset,
            vec![
                atom_ligand(2),
                virtual_ligand(0),
                atom_ligand(3),
                virtual_ligand(1),
            ],
        );

        assert_eq!(
            projected_cis_trans_bonds(&molecule),
            [CisTransBond {
                bond: 0,
                first_ligand: 2,
                second_ligand: 3,
                relation,
            }]
        );
    }

    #[rstest]
    #[case::virtual_first(
        vec![virtual_ligand(0), atom_ligand(2), atom_ligand(3), virtual_ligand(1)],
        SideRelation::OppositeSide
    )]
    #[case::reversed_endpoints(
        vec![atom_ligand(3), virtual_ligand(1), atom_ligand(2), virtual_ligand(0)],
        SideRelation::SameSide
    )]
    fn test_cis_trans_bond_reframing(
        #[case] ligands: Vec<StereoLigand>,
        #[case] relation: SideRelation,
    ) {
        let molecule = stereo_molecule(StereoKind::CisTrans, StereoCoset::Lit(0), ligands);

        assert_eq!(
            projected_cis_trans_bonds(&molecule),
            [CisTransBond {
                bond: 0,
                first_ligand: 2,
                second_ligand: 3,
                relation,
            }]
        );
    }

    #[rstest]
    #[case::undetermined(StereoKind::CisTrans, StereoCoset::Undetermined)]
    #[case::set(
        StereoKind::CisTrans,
        StereoCoset::lit_set([0, 1])
    )]
    #[case::term(
        StereoKind::CisTrans,
        StereoCoset::term(StereoTerm::var("configuration"))
    )]
    #[case::axial(StereoKind::Axial, StereoCoset::Lit(0))]
    fn test_cis_trans_bond_omission(#[case] kind: StereoKind, #[case] coset: StereoCoset) {
        let molecule = stereo_molecule(
            kind,
            coset,
            vec![
                atom_ligand(2),
                virtual_ligand(0),
                atom_ligand(3),
                virtual_ligand(1),
            ],
        );

        assert_eq!(projected_cis_trans_bonds(&molecule), []);
    }

    #[rstest]
    fn test_cis_trans_bond_fully_undetermined() {
        let molecule = stereo_molecule_with_attributes(
            StereoBondForm::default(),
            vec![
                atom_ligand(2),
                virtual_ligand(0),
                atom_ligand(3),
                virtual_ligand(1),
            ],
        );

        assert_eq!(projected_cis_trans_bonds(&molecule), []);
    }

    fn projected_cis_trans_bonds(molecule: &Molecule) -> Vec<CisTransBond> {
        molecule
            .stereo_bonds()
            .iter()
            .filter_map(|stereo| cis_trans_bond(molecule, stereo))
            .collect()
    }

    fn stereo_molecule(
        kind: StereoKind,
        coset: StereoCoset,
        ligands: Vec<StereoLigand>,
    ) -> Molecule {
        stereo_molecule_with_attributes(StereoBondForm::new(kind, coset), ligands)
    }

    fn stereo_molecule_with_attributes(
        attributes: StereoBondForm,
        ligands: Vec<StereoLigand>,
    ) -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: [Element::C, Element::C, Element::F, Element::Cl]
                .into_iter()
                .map(AtomForm::from_element)
                .collect(),
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(2)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(1), AtomId(3), BondForm::from_order(1)),
            ],
            stereo_bonds: vec![(BondId(0), ligands, attributes)],
            ..Default::default()
        })
    }

    fn atom_ligand(atom: u32) -> StereoLigand {
        StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)
    }

    fn virtual_ligand(atom: u32) -> StereoLigand {
        StereoLigand::new(AtomId(atom), StereoLigandKind::ImplicitHydrogen)
    }
}
