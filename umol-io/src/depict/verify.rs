//! Verification of a supplied molecule layout against the molecule it depicts.

use umol_geometric_core::{
    finite_difference, normalizable_direction, same_side_of_axis, Point2D, Point3D,
};
use umol_graph_ir::ir::{AtomId, CisTransConfiguration, Entity, Molecule};

use super::molecule::{tetrahedral_wedges, MoleculeDepictionError};
use crate::layout::stereo::{cis_trans_site, CisTransSite};
use crate::layout::MoleculeLayout;

/// Checks frame agreement, definite cis/trans agreement, finite derived geometry, and tetrahedral
/// wedge selection, in order.
pub(crate) fn verify_molecule_layout(
    molecule: &Molecule,
    layout: &MoleculeLayout,
) -> Result<(), MoleculeDepictionError> {
    layout.check_frame(molecule)?;
    for stereo in molecule.stereo_bonds().iter() {
        if let Some(site) = cis_trans_site(molecule, stereo) {
            check_cis_trans(layout, site)?;
        }
    }
    check_finite_geometry(molecule, layout)?;
    tetrahedral_wedges(molecule, layout).map(drop)
}

fn check_cis_trans(
    layout: &MoleculeLayout,
    site: CisTransSite,
) -> Result<(), MoleculeDepictionError> {
    let [site_0, site_1] = site.site.map(|atom| point3(position(layout, atom)));
    let ligands = site.ligands.map(|atom| point3(position(layout, atom)));
    let found = match same_side_of_axis(site_0, site_1, ligands[0], ligands[1]) {
        Some(true) => CisTransConfiguration::Z,
        Some(false) => CisTransConfiguration::E,
        None => {
            return Err(MoleculeDepictionError::CisTransDegenerate {
                bond: site.bond,
                ligand: degenerate_ligand(site_0, site_1, site.ligands, ligands),
            })
        }
    };
    if found != site.configuration {
        return Err(MoleculeDepictionError::CisTransMismatch {
            bond: site.bond,
            expected: site.configuration,
            found,
        });
    }
    Ok(())
}

fn degenerate_ligand(
    site_0: Point3D,
    site_1: Point3D,
    ligands: [AtomId; 2],
    points: [Point3D; 2],
) -> AtomId {
    ligands
        .into_iter()
        .zip(points)
        .find(|&(_, point)| same_side_of_axis(site_0, site_1, point, point).is_none())
        .map_or(ligands[0], |(ligand, _)| ligand)
}

fn check_finite_geometry(
    molecule: &Molecule,
    layout: &MoleculeLayout,
) -> Result<(), MoleculeDepictionError> {
    for bond in molecule.bonds().iter() {
        let [first, second] = bond.atom_ids().map(|atom| position(layout, atom));
        if normalizable_direction(first, second).is_none() {
            return Err(MoleculeDepictionError::NonFiniteGeometry {
                entity: Entity::Bond(bond.id),
            });
        }
    }
    let Some(origin) = layout
        .positions()
        .iter()
        .copied()
        .reduce(|min, position| Point2D::new(min.x.min(position.x), min.y.min(position.y)))
    else {
        return Ok(());
    };
    for atom in molecule.atoms().iter() {
        if finite_difference(origin, position(layout, atom.id)).is_none() {
            return Err(MoleculeDepictionError::NonFiniteGeometry {
                entity: Entity::Atom(atom.id),
            });
        }
    }
    Ok(())
}

fn position(layout: &MoleculeLayout, atom: AtomId) -> Point2D {
    *layout
        .position(atom)
        .expect("frame agreement establishes every graph-IR atom position")
}

fn point3(point: Point2D) -> Point3D {
    Point3D::new(point.x, point.y, 0.0)
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{BondId, StereoAtomId};
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::depict::Depict;
    use crate::layout::MoleculeLayoutError;

    const TRANS_BUTENE: &str = r#"{:atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#;
    const CIS_BUTENE: &str = r#"{:atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
        :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]}"#;

    fn layout(positions: &[[f64; 2]]) -> MoleculeLayout {
        MoleculeLayout::try_new(positions.iter().map(|&[x, y]| Point2D::new(x, y)).collect())
            .unwrap()
    }

    #[rstest]
    #[case::trans(TRANS_BUTENE)]
    #[case::cis(CIS_BUTENE)]
    fn test_verify_generated_layout(#[case] input: &str) {
        let molecule = mol_dsl!(input);
        let generated = molecule.layout().unwrap();

        assert_eq!(molecule.verify_layout(&generated), Ok(()));
    }

    #[rstest]
    #[case::trans_drawn_as_cis(
        TRANS_BUTENE,
        [[0.0, 1.0], [1.0, 0.0], [2.0, 0.0], [3.0, 1.0]],
        CisTransConfiguration::E,
        CisTransConfiguration::Z
    )]
    #[case::cis_drawn_as_trans(
        CIS_BUTENE,
        [[0.0, 1.0], [1.0, 0.0], [2.0, 0.0], [3.0, -1.0]],
        CisTransConfiguration::Z,
        CisTransConfiguration::E
    )]
    fn test_verify_cis_trans_mismatch(
        #[case] input: &str,
        #[case] positions: [[f64; 2]; 4],
        #[case] expected: CisTransConfiguration,
        #[case] found: CisTransConfiguration,
    ) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            molecule.verify_layout(&layout(&positions)),
            Err(MoleculeDepictionError::CisTransMismatch {
                bond: BondId(1),
                expected,
                found,
            })
        );
    }

    #[rstest]
    #[case::trans_matches(TRANS_BUTENE, [[0.0, 1.0], [1.0, 0.0], [2.0, 0.0], [3.0, -1.0]])]
    #[case::cis_matches(CIS_BUTENE, [[0.0, 1.0], [1.0, 0.0], [2.0, 0.0], [3.0, 1.0]])]
    fn test_verify_cis_trans_agreement(#[case] input: &str, #[case] positions: [[f64; 2]; 4]) {
        let molecule = mol_dsl!(input);

        assert_eq!(molecule.verify_layout(&layout(&positions)), Ok(()));
    }

    #[rstest]
    #[case::first_ligand_on_axis([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 1.0]], AtomId(0))]
    #[case::second_ligand_on_axis([[0.0, 1.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]], AtomId(3))]
    #[case::zero_length_site_bond([[0.0, 1.0], [1.0, 0.0], [1.0, 0.0], [3.0, -1.0]], AtomId(0))]
    fn test_verify_cis_trans_degenerate(#[case] positions: [[f64; 2]; 4], #[case] ligand: AtomId) {
        let molecule = mol_dsl!(TRANS_BUTENE);

        assert_eq!(
            molecule.verify_layout(&layout(&positions)),
            Err(MoleculeDepictionError::CisTransDegenerate {
                bond: BondId(1),
                ligand,
            })
        );
    }

    #[rstest]
    #[case::coincident_bonded_atoms(
        r#"{:atoms ["C" "O" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#,
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 0.0]],
        Entity::Bond(BondId(1))
    )]
    #[case::overflowing_bond_difference(
        r#"{:atoms ["C" "O"] :bonds [[0 1 "1"]]}"#,
        vec![[1e308, 0.0], [-1e308, 0.0]],
        Entity::Bond(BondId(0))
    )]
    #[case::nonliteral_bond_order(
        r#"{:atoms ["C" "O"] :bonds [[0 1 "*"]]}"#,
        vec![[0.5, 0.5], [0.5, 0.5]],
        Entity::Bond(BondId(0))
    )]
    #[case::overflowing_extent(
        r#"{:atoms ["C" "O"] :bonds []}"#,
        vec![[-1e308, 0.0], [1e308, 0.0]],
        Entity::Atom(AtomId(1))
    )]
    fn test_verify_non_finite_geometry(
        #[case] input: &str,
        #[case] positions: Vec<[f64; 2]>,
        #[case] entity: Entity,
    ) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            molecule.verify_layout(&layout(&positions)),
            Err(MoleculeDepictionError::NonFiniteGeometry { entity })
        );
    }

    #[rstest]
    fn test_verify_layout_frame() {
        let molecule = mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#);

        assert_eq!(
            molecule.verify_layout(&layout(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])),
            Err(MoleculeDepictionError::LayoutFrame(
                MoleculeLayoutError::FrameSizeMismatch {
                    molecule_atom_count: 2,
                    layout_atom_count: 3,
                }
            ))
        );
    }

    #[rstest]
    fn test_verify_tetrahedral_geometry() {
        let molecule = mol_dsl!(
            r#"{:atoms ["C" "F" "Cl" "Br" "I"]
                :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#
        );
        let collinear = layout(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0], [4.0, 0.0]]);

        assert_eq!(
            molecule.verify_layout(&collinear),
            Err(MoleculeDepictionError::TetrahedralGeometry {
                stereo_atom: StereoAtomId(0),
            })
        );
    }
}
